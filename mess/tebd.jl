using OMEinsum
using Yao
using LightGraphs
using Yao.ConstGate: P1
using LinearAlgebra: svd
using LightGraphs: SimpleEdge

function cterm(C::Real)
    return Matrix(C*kron(P1, P1))
end

function zterm(Δ::Real)
    return Matrix(Δ*Z)
end

function xterm(Ω::Real)
    return Matrix(Ω*X)
end

struct PEPS{T, LT<:Union{Int,Char}}
    physical_labels::Vector{LT}
    virtual_labels::Vector{LT}

    vertex_labels::Vector{Vector{LT}}
    vertex_tensors::Vector{<:AbstractArray{T}}
    # bond labels are same as virtual_labels
    bond_tensors::Vector{<:AbstractVector{T}}
    max_index::LT
end

function state(peps::PEPS)
    code = EinCode((Tuple.(peps.vertex_labels)..., Tuple.(peps.virtual_labels)...), Tuple(peps.physical_labels))
    size_dict = OMEinsum.get_size_dict(peps.vertex_labels, peps.vertex_tensors)
    optcode = optimize_greedy(code, size_dict)
    optcode(peps.vertex_tensors..., peps.bond_tensors...)
end
getvlabel(peps::PEPS, i::Int) = peps.vertex_labels[i]
getblabel(peps::PEPS, i::Int) = [peps.virtual_labels[i]]
getphysicallabel(peps::PEPS, i::Int) = peps.physical_labels[i]
newlabel(peps, offset) = peps.max_index + offset
findbondtensor(peps, b) = (peps.bond_tensors[findall(==(b), peps.virtual_labels)[]])

function applyA!(peps, g, Δ::Real, Ω::Real, δt::Real)
    for v in vertices(g)
        hi = xterm(Ω)*zterm(Δ)
        te = exp(-im*δt*hi)
        apply_onsite!(peps, v, te)
    end
end

function applyB!(peps, g, C::Real)
    for (i,j) in edges(g)
        hi = cterm(C)
        te = reshape(exp(-im*δt*hi), 2, 2, 2, 2)
        apply_onbond!(peps, i, j, te)
    end
end

function apply_onsite!(peps::PEPS{T,LT}, i, mat::AbstractMatrix) where {T,LT}
    @assert size(mat, 1) == size(mat, 2)
    ti = peps.vertex_tensors[i]
    old = (getvlabel(peps, i)...,)
    mlabel = (getphysicallabel(peps, i), newlabel(peps, 1))
    peps.vertex_tensors[i] = EinCode((old, mlabel), replace(old, mlabel[1]=>mlabel[2]))(ti, mat)
    return peps
end

function apply_onbond!(peps::PEPS, i, j, mat::AbstractArray{T,4}, D) where T
    ti, tj = peps.vertex_tensors[i], peps.vertex_tensors[j]
    li, lj = getvlabel(peps, i), getvlabel(peps, j)
    shared_label = li ∩ lj; @assert length(shared_label) == 1
    only_left, only_right = setdiff(li, lj), setdiff(lj, li)
    lij = ((only_left ∪ only_right)...,)

    # multiple S tensors
    for b in li
        b == getphysicallabel(peps, i) && continue
        ti = _apply_vec(peps, ti, li, findbondtensor(peps, b), b)
    end
    for b in only_right
        b == getphysicallabel(peps, j) && continue
        tj = _apply_vec(peps, tj, lj, findbondtensor(peps, b), b)
    end
    tij = EinCode(((li...,), (lj...,)), (lij...,))(ti, tj)
    lijkl = (getphysicallabel(peps, i), getphysicallabel(peps, j), newlabel(peps, 1), newlabel(peps, 2))
    lkl = replace(lij, lijkl[1]=>lijkl[3], lijkl[2]=>lijkl[4])
    tkl = EinCode(((lij...,), lijkl), (lkl...,))(tij, mat)

    # SVD and truncate
    sl = [size(ti, findall(==(l), li)[]) for l in only_left]
    sr = [size(tj, findall(==(l), lj)[]) for l in only_right]
    U, S, Vt = svd(reshape(tkl, prod(sl), prod(sr)))
    if D !== nothing
        S = S[1:D]
        U = U[:,1:D]
        Vt = conj.(Vt[:,1:D])
        println("truncation error is $(sum(S[D+1:end]))")
    end
    
    # reshape back
    ti_ = EinCode(((only_left..., shared_label[]),), (li...,))(reshape(U, (sl..., size(U,2))))
    tj_ = EinCode(((only_right..., shared_label[]),), (lj...,))(reshape(Vt, (sr..., size(Vt,2))))

    # devide S tensors
    for b in li
        b == getphysicallabel(peps, i) && continue
        ti_ = _apply_vec(peps, ti_, li, 1 ./ replace(findbondtensor(peps, b), T(0)=>T(1e-20)), b)  # assume S being positive
    end
    for b in only_right
        b == getphysicallabel(peps, j) && continue
        tj_ = _apply_vec(peps, tj_, lj, 1 ./ replace(findbondtensor(peps, b), T(0)=>T(1e-20)), b)
    end

    # update tensors
    peps.vertex_tensors[i] = ti_
    peps.vertex_tensors[j] = tj_
    peps.bond_tensors[findall(==(shared_label[]), peps.virtual_labels)[]] = S
    return peps
end

function _apply_vec(peps, t, l, v, b)
    return EinCode(((l...,), (b,)), (l...,))(t, v)
end

function PEPS(vertex_labels::AbstractVector{<:AbstractVector{LT}}, vertex_tensors::Vector{<:AbstractArray{T}},
        virtual_labels::AbstractVector{LT}, bond_tensors::Vector{<:AbstractVector}) where {LT,T}
    physical_labels = [vl[findall(∉(virtual_labels), vl)[]] for vl in vertex_labels]
    max_ind = max(maximum(physical_labels), maximum(virtual_labels))
    PEPS(physical_labels, virtual_labels, vertex_labels, vertex_tensors, bond_tensors, max_ind)
end

function peps_zero_state(::Type{T}, g::SimpleGraph, D::Int) where T
    virtual_labels = collect(nv(g)+1:nv(g)+ne(g))
    vertex_labels = Vector{Int}[]
    vertex_tensors = Array{T}[]
    edge_map = Dict(zip(edges(g), virtual_labels))
    for i=1:nv(g)
        push!(vertex_labels, [i,[get(edge_map, SimpleEdge(i,nb), get(edge_map,SimpleEdge(nb,i),0)) for nb in neighbors(g, i)]...])
        t = zeros(T, 2, fill(D, degree(g, i))...)
        t[1] = 1
        push!(vertex_tensors, t)
    end
    if any(vl->any(iszero, vl), vertex_labels)
        error("incorrect input labels1")
    end
    bond_tensors = [ones(T, D) for _=1:ne(g)]
    return PEPS(vertex_labels, vertex_tensors, virtual_labels, bond_tensors)
end
using Test

@testset "initial state" begin
    g = SimpleGraph(5)
    for (i,j) in [(1,2), (1,3), (2,4), (2,5), (3,4), (3,5)]
        add_edge!(g, i, j)
    end
    peps = peps_zero_state(ComplexF64, g, 2)
    @test vec(state(peps)) ≈ (x=zeros(ComplexF64, 1<<5); x[1]=1; x)
    @test newlabel(peps, 2) == 11+2
    @test newlabel(peps, 4) == 11+4
    @test getphysicallabel(peps, 2) == 2
    @test length(getvlabel(peps, 2)) == 4
    @test getblabel(peps, 2) == [7]
    apply_onsite!(peps, 1, [0 1; 1 0])
    @test vec(state(peps)) ≈ (x=zeros(ComplexF64, 1<<5); x[2]=1; x)
    apply_onbond!(peps, 1, 2, reshape(Matrix(ConstGate.CNOT), 2, 2, 2, 2), 2)
    @test vec(state(peps)) ≈ (x=zeros(ComplexF64, 1<<5); x[4]=1; x)
end