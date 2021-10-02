using OMEinsum
using Yao
using LightGraphs
using Yao.ConstGate: P1
using LinearAlgebra: svd

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
    #graph::SimpleGraph{Int}
    code::EinCode
    physical_labels::Vector{LT}
    tensors::Vector{AbstractArray{T}}
    max_index::LT
end

state(peps::PEPS) = peps.code(peps.tensors...)
getlabel(peps::PEPS, i::Int) = getixs(peps.code)[i]
getphysicallabel(peps::PEPS, i::Int) = peps.physical_labels[i]
newlabel(peps, offset) = peps.max_index + offset

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
    ti = peps.tensors[i]
    old = getlabel(peps, i)
    mlabel = (getphysicallabel(peps, i), newlabel(peps, 1))
    peps.tensors[i] = EinCode((old, mlabel), replace(old, mlabel[1]=>mlabel[2]))(ti, mat)
    return peps
end

function apply_onbond!(peps::PEPS, i, j, mat::AbstractArray{T,4}) where T
    ti, tj = peps.tensors[i], peps.tensors[j]
    li, lj = getlabel(peps, i), getlabel(peps, j)
    shared_label = li ∩ lj; @assert length(shared_label) == 1
    only_left, only_right = setdiff(li, lj), setdiff(lj, li)
    lij = ((only_left ∪ only_right)...,)
    tij = EinCode((li, lj), lij)(ti, tj)
    lijkl = (getphysicallabel(peps, i), getphysicallabel(peps, j), newlabel(peps, 1), newlabel(peps, 2))
    lkl = replace(lij, lijkl[1]=>lijkl[3], lijkl[2]=>lijkl[4])
    tkl = EinCode((lij, lijkl), lkl)(tij, mat)
    # SVD and truncate
    U, S, V = svd(reshape(tkl, prod(sl), prod(sr)))
    sqrtS = Diagonal(sqrt(S[1:D]))
    tl = U[:,1:D] * sqrtS
    tr = sqrtS * V[1:D,:]
    println("truncation error is $(sum(S[D+1:end]))")

    # reshape back
    peps.tensors[i] = EinCode(((only_left..., shared_label[1]),), li)(reshape(tl, slnew))
    peps.tensors[j] = EinCode(((shared_label[1], only_right...),), lj)(reshape(tr, srnew))
    #L, R = svd_compress(tkl, lkl, left, right, D)
    return peps
end