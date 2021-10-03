using LinearAlgebra: I, eigvals, normalize
using Random
n1, n2 = 1000, 20
ρ = 0.01
E = 0.0
h11 = E .* Matrix(I, n1, n1)
h22 = zeros(n2, n2)
Random.seed!(2)
h12 = (rand(n1, n2) .< ρ)
h = [h11 h12; h12' h22]

evals = eigvals(h)
gap = evals[n2+1] - evals[n2]
estimate = 2*normalize(ones(n1))' * h12 * normalize(ones(n2))
@show gap estimate