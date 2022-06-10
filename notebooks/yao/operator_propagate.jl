### A Pluto.jl notebook ###
# v0.19.8

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 69a8bfdc-e834-11ec-2a4d-b31fc1acc575
using Pkg; Pkg.develop("Yao")

# ╔═╡ 5b498022-64d6-42f1-8dde-e89a3941d76a
using Yao, YaoPlots

# ╔═╡ 1fd55c3a-cb44-4202-b7fd-e98f311c7fba
using GenericTensorNetworks, Graphs

# ╔═╡ 8752afba-23b7-4c89-add9-b3e80dba29aa
using PlutoUI

# ╔═╡ 93a3a592-338c-45fa-af04-9ed86e205043
let
	gadget = @bind Δ Slider(0.0:0.01:3.0, default=0.4, show_value=true)
	md"Δ = $gadget"
end

# ╔═╡ 5c33f5e2-4e3a-4714-8ed2-7537c7bbbdb7
let
	gadget = @bind Ω Slider(0.0:0.01:3.0, default=0.4, show_value=true)
	md"Ω = $gadget"
end

# ╔═╡ d3127510-7b7f-4108-b679-44d8bea640bc
graph = random_diagonal_coupled_graph(8, 8, 0.8)

# ╔═╡ 0c4c03f2-76b2-4c13-b551-1aae3a5b4c87
xterms = Ω * sum([put(nv(graph), i=>X) for i=1:nv(graph)]);

# ╔═╡ ee47c7d3-06ad-4cbd-92de-95aad1f0ee2f
zterms = Δ * sum([put(nv(graph), i=>Z) for i=1:nv(graph)]);

# ╔═╡ 4d643785-ece2-4ecc-a292-d6937930f9af
function projector(::Type{T}, configs::AbstractVector{<:DitStr{D,N}}) where {T,D,N}
	mask = zeros(T, D^N)
	for c in configs
		mask[Int(c)+1] = one(T)
	end
	return matblock(Diagonal(mask); tag="Projector")
end

# ╔═╡ 9c8f1c0a-9ae6-472d-aa75-702d6ad61c32
function pxp(graph, i::Int)
	n = nv(graph)
	P = chain([put(n, (i,j)=>projector(ComplexF64, [bit"10", bit"01", bit"00"])) for j in neighbors(graph, i)])
	P * put(n, i=>X) * P
end

# ╔═╡ b43d5fe3-1190-45f3-8239-e16893f0100f
pxpterms = Ω * sum([pxp(graph, i) for i=1:nv(graph)]);

# ╔═╡ 999ad599-fb65-4afa-8a3d-cc471b61cb4f
projector(ComplexF64, [bit"10", bit"01", bit"00"])

# ╔═╡ 4d499494-82dd-453f-9633-798565773048
pxp(graph, 1)

# ╔═╡ 3f7f60cb-5e44-43ec-b80b-b3ea00e5663d
pxpterms[:,zero(BitStr{nv(graph),Int})]

# ╔═╡ 5141f864-1b54-4514-b2d0-eb1cd56f920b
pxpterms[:,zero(BitStr{nv(graph),Int128})]

# ╔═╡ d8da9fc3-9715-450a-8d34-752c7389c87b
table = pxpterms[:,pxpterms[:,zero(BitStr{nv(graph),Int128})]]

# ╔═╡ f5c50fc9-d3d3-4ebe-a739-a49f11c5a397
cleanup(table; zero_threshold=0.33)

# ╔═╡ d9afea06-2d11-475d-925b-44ccf96b551a
md"## Suboptimal states"

# ╔═╡ 3cafacc2-fd0b-412c-b4ff-7836c840ce5f
optimal = solve(IndependentSet(graph), SingleConfigMax())[].c.data

# ╔═╡ f185e12c-44bf-4ce8-a2f6-8c32783aab1a
DitStr{2}(optimal)

# ╔═╡ f9bc40b9-e247-4f55-ae27-353a838aeef8
totalsize = solve(IndependentSet(graph), CountingAll())[]

# ╔═╡ 79e37614-7dd6-4950-9ed9-7f97457ecdef
@time let
	j = EntryTable([DitStr{2}(optimal)], [1.0+0im])
	for i=1:5
		@show i
		j = pxpterms[:, j]
	end
	length(j)
end

# ╔═╡ 0221e1c8-5a46-40be-ab07-2159215bd9dd
pxpterms.content |> length

# ╔═╡ 86034570-e6c6-4a4b-a472-6f690ae1fc7c
md"## Speed up!"

# ╔═╡ 242f1a78-fcef-4870-8def-28d5915cf34b
flt = Optimise.simplify(pxpterms; rules=[Optimise.eliminate_nested])

# ╔═╡ 128716f6-0b7e-44ba-aa41-124e983818c3
function squeeze_locations(pb::PutBlock, locs)
	put(length(locs), map(loc->findfirst(==(loc), locs), pb.locs)=>pb.content)
end

# ╔═╡ aa83395c-feec-4ae0-997e-b624c12b2ff7
begin
	function rule_chain2put(c::ChainBlock)
		locs = occupied_locs(c)
		put(nqudits(c), locs=>matblock(mat(chain(squeeze_locations.(subblocks(c), Ref(locs)))); tag="PXP"))
	end
	rule_chain2put(c::AbstractBlock) = chsubblocks(c, rule_chain2put.(subblocks(c)))
end

# ╔═╡ 028508b5-92f9-488d-a549-293a054dd6d9
rule_chain2put(flt.content[1])

# ╔═╡ 0ccd384a-d383-45c2-876b-f9bc528946ba
smpl = rule_chain2put(flt)

# ╔═╡ 3031a8e8-51b4-4e56-bcdd-c9a86d25e1a5
@time let
	j = EntryTable([DitStr{2}(optimal)], [1.0+0im])
	for i=1:5
		@show i
		j = smpl[:, j]
	end
	length(j)
end

# ╔═╡ Cell order:
# ╠═69a8bfdc-e834-11ec-2a4d-b31fc1acc575
# ╠═5b498022-64d6-42f1-8dde-e89a3941d76a
# ╠═1fd55c3a-cb44-4202-b7fd-e98f311c7fba
# ╠═8752afba-23b7-4c89-add9-b3e80dba29aa
# ╟─93a3a592-338c-45fa-af04-9ed86e205043
# ╟─5c33f5e2-4e3a-4714-8ed2-7537c7bbbdb7
# ╠═d3127510-7b7f-4108-b679-44d8bea640bc
# ╠═0c4c03f2-76b2-4c13-b551-1aae3a5b4c87
# ╠═ee47c7d3-06ad-4cbd-92de-95aad1f0ee2f
# ╠═4d643785-ece2-4ecc-a292-d6937930f9af
# ╠═9c8f1c0a-9ae6-472d-aa75-702d6ad61c32
# ╠═b43d5fe3-1190-45f3-8239-e16893f0100f
# ╠═999ad599-fb65-4afa-8a3d-cc471b61cb4f
# ╠═4d499494-82dd-453f-9633-798565773048
# ╠═3f7f60cb-5e44-43ec-b80b-b3ea00e5663d
# ╠═5141f864-1b54-4514-b2d0-eb1cd56f920b
# ╠═d8da9fc3-9715-450a-8d34-752c7389c87b
# ╠═f5c50fc9-d3d3-4ebe-a739-a49f11c5a397
# ╟─d9afea06-2d11-475d-925b-44ccf96b551a
# ╠═3cafacc2-fd0b-412c-b4ff-7836c840ce5f
# ╠═f185e12c-44bf-4ce8-a2f6-8c32783aab1a
# ╠═f9bc40b9-e247-4f55-ae27-353a838aeef8
# ╠═79e37614-7dd6-4950-9ed9-7f97457ecdef
# ╠═0221e1c8-5a46-40be-ab07-2159215bd9dd
# ╟─86034570-e6c6-4a4b-a472-6f690ae1fc7c
# ╠═242f1a78-fcef-4870-8def-28d5915cf34b
# ╠═aa83395c-feec-4ae0-997e-b624c12b2ff7
# ╠═128716f6-0b7e-44ba-aa41-124e983818c3
# ╠═028508b5-92f9-488d-a549-293a054dd6d9
# ╠═0ccd384a-d383-45c2-876b-f9bc528946ba
# ╠═3031a8e8-51b4-4e56-bcdd-c9a86d25e1a5
