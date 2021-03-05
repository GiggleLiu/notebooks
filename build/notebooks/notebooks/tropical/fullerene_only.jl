### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 63bae6c0-7a00-11eb-1517-bdf4ecf20377
using SimpleTensorNetworks, TropicalNumbers, Random

# ╔═╡ c547d8c4-7a00-11eb-185a-353cac3a747b
md"## The Song Shan Lake Spring School (SSSS) Challendge"

# ╔═╡ c7953d3a-7a00-11eb-1338-5d69526ea645
md"""
[You may find more than 10 solutions in our github repo](https://github.com/QuantumBFS/SSSS/blob/master/Challenge.md)
"""

# ╔═╡ 4fa6e6a2-7a00-11eb-08c8-09659b05f157
function ising_bondtensor(::Type{T}, J) where T
	e = T(J)
	e_ = T(-J)
	[e e_; e_ e]
end

# ╔═╡ 5fb276c4-7a00-11eb-3c90-390e289777cd
function ising_vertextensor(::Type{T}, n::Int, h) where T
	res = zeros(T, fill(2, n)...)
	res[1] = T(h)
	res[end] = T(-h)
	return res
end

# ╔═╡ 6ecd8b12-7a00-11eb-0330-9b591bcfa756
edges = [
	 	1=>10, 1=>41, 1=>59, 2=>12, 2=>42, 2=>60, 3=>6, 3=>
        43, 3=>57, 4=>8, 4=>44, 4=>58, 5=>13, 5=>56, 5=>
        57, 6=>10, 6=>31, 7=>14, 7=>56, 7=>58, 8=>12, 8=>
        32, 9=>23, 9=>53, 9=>59, 10=>15, 11=>24, 11=>53, 11=>
        60, 12=>16, 13=>14, 13=>25, 14=>26, 15=>27, 15=>
        49, 16=>28, 16=>50, 17=>18, 17=>19, 17=>54, 18=>
        20, 18=>55, 19=>23, 19=>41, 20=>24, 20=>42, 21=>
        31, 21=>33, 21=>57, 22=>32, 22=>34, 22=>58, 23=>
        24, 25=>35, 25=>43, 26=>36, 26=>44, 27=>51, 27=>
        59, 28=>52, 28=>60, 29=>33, 29=>34, 29=>56, 30=>
        51, 30=>52, 30=>53, 31=>47, 32=>48, 33=>45, 34=>
        46, 35=>36, 35=>37, 36=>38, 37=>39, 37=>49, 38=>
        40, 38=>50, 39=>40, 39=>51, 40=>52, 41=>47, 42=>
        48, 43=>49, 44=>50, 45=>46, 45=>54, 46=>55, 47=>
        54, 48=>55
];

# ╔═╡ bb691266-7a00-11eb-2da9-f1c63699d905
md"construct tensor network by assigning labels `(edge index, boolean)` to tensors, where the boolean identifies whether this label correspond to the source node or the destination side of the edge. The resulting tensor network contains 90 edge tensors and 60 vertex tensors."

# ╔═╡ 7c8564e6-7a00-11eb-2667-3b8fbf1a017d
tn = TensorNetwork(vcat(
	[LabeledTensor(ising_vertextensor(CountingTropical{Float64}, 3, 0.0), [(i,j==e.first) for (i,e) in enumerate(edges) if j ∈ e]) for j=1:60],  # vertex
	[LabeledTensor(ising_bondtensor(CountingTropical{Float64}, -1.0), [(i,true),(i,false)]) for i = 1:length(edges)]
));

# ╔═╡ 9144964a-7a00-11eb-3493-47c93da5cb19
# find a good contraction order by greedy search
tcs, scs, order = (Random.seed!(2); trees_greedy(tn; strategy="min_reduce"));

# ╔═╡ a2baab12-7a00-11eb-1d12-65de2c78cca3
md"time complexity = $(round(log2sumexp2(tcs); sigdigits=4)), space complexity = $(round(maximum(scs); sigdigits=4))"

# ╔═╡ 99303312-7a00-11eb-2d20-57acdf9ad2a8
SimpleTensorNetworks.contract(tn, order[]).array[]

# ╔═╡ Cell order:
# ╠═63bae6c0-7a00-11eb-1517-bdf4ecf20377
# ╟─c547d8c4-7a00-11eb-185a-353cac3a747b
# ╟─c7953d3a-7a00-11eb-1338-5d69526ea645
# ╠═4fa6e6a2-7a00-11eb-08c8-09659b05f157
# ╠═5fb276c4-7a00-11eb-3c90-390e289777cd
# ╠═6ecd8b12-7a00-11eb-0330-9b591bcfa756
# ╟─bb691266-7a00-11eb-2da9-f1c63699d905
# ╠═7c8564e6-7a00-11eb-2667-3b8fbf1a017d
# ╠═9144964a-7a00-11eb-3493-47c93da5cb19
# ╟─a2baab12-7a00-11eb-1d12-65de2c78cca3
# ╠═99303312-7a00-11eb-2d20-57acdf9ad2a8
