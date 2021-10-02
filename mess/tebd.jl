using OMEinsum
using Yao
using LightGraphs
using Yao.ConstGate: P1

function cterm(C::Real)
    return Matrix(C*kron(P1, P1))
end

function zterm(Δ::Real)
    return Matrix(Δ*Z)
end

function xterm(Ω::Real)
    return Matrix(Ω*X)
end

struct PEPS{T}
    #graph::SimpleGraph{Int}
    code
    tensors::Vector{AbstractArray{T}}
end

state(peps::PEPS) = peps.code(peps.tensors...)

getlabel(peps::PEPS, i::Int) = getixs(code)[i]

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

function apply_onsite!(peps::PEPS, i, label, mat::AbstractMatrix)
    @assert size(mat, 1) == size(mat, 2)
    ti = peps.tensors[i]
    m = ndims(ti)
    lbs = getlabel(peps, i)
    ilabel = findfirst(==(label), lbs)
    if ilabel === nothing
        error("label $label does not exist in tensor $lbs")
    else
        old = collect(1:m)
        new = copy(old)
        new[ilabel] = m+1
        ti .= EinCode(((old...,), (ilabel, m+1)), (new...,))(ti, mat)
        return peps
    end
end

function apply_onbond!(peps::PEPS, i, j, mat::AbstractArray{T,4}) where T
    ti, tj = peps.tensors[i], peps.tensors[j]
    li, lj = getlabel(peps, i), getlabel(peps, j)
    tij = EinCode((li, lj), lij)(ti, tj)
    tkl = EinCode((lij, lijkl), lkl)(tij, mat)
    L, R = svd_compress(tkl, lkl, left, right)
    ti .= EinCode(ntuple(identity, m), (m, m+1))(ti, mat)
end

function svd_compress(t, label, left, right, D::Int)
    # permute input tensors
    lefto = left ∩ t
    righto = right ∩ t
    ly = (lefto ∪ righto)
    y = EinCode((label,), (ly...,))(t)

    # svd and truncate
    U, S, V = svd(reshape(y, Dl, Dr))
    sqrtS = Diagonal(sqrt(S[1:D]))
    tl = U[:,1:D] * sqrtS
    tr = sqrtS * V[1:D,:]

    # reshape back
    left_ = reshape(tl, sl)
    right_ = reshape(tr, sr)
    println("truncation error is $(sum(S[D+1:end]))")
    return EinCode((left_,), left)(reshape(tl,sl_)), EinCode((right_,), right)(tr,sr_)
end