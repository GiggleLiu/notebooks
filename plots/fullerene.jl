### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 63bae6c0-7a00-11eb-1517-bdf4ecf20377
using Revise, SimpleTensorNetworks, TropicalNumbers, Random

# ╔═╡ 4dac28aa-c273-454f-981e-0087e3dbd939
using LightGraphs, Viznet, Compose, PlutoUI

# ╔═╡ a2edc836-42fd-4723-822b-540019980726
using ForwardDiff

# ╔═╡ 939d981f-5364-46e5-90a6-cb23d7cf8f7f
using Rotations

# ╔═╡ e4ff4e4e-6584-4a7b-839c-892647365201
using StaticArrays, CoordinateTransformations

# ╔═╡ 208080a7-392d-43af-aebd-3f1d8e5af325
using GraphPlot

# ╔═╡ 5953b211-2f32-4e1c-9ebc-77ce07c5bcc5
using Statistics

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

# ╔═╡ b65e6868-efdc-4b34-938a-6d5fbbad1e16
function solve(edges, hs::AbstractVector{T}) where T
	tn = TensorNetwork(vcat(
		[LabeledTensor(ising_vertextensor(Tropical{T}, 3, hs[j]), [(i,j==e.first) for (i,e) in enumerate(edges) if j ∈ e]) for j=1:60],  # vertex
		[LabeledTensor(ising_bondtensor(Tropical{T}, T(-1.0)), [(i,true),(i,false)]) for i = 1:length(edges)]
	))
	tcs, scs, order = trees_greedy(tn; strategy="min_reduce")
	SimpleTensorNetworks.contract(tn, order[]).array[]
end

# ╔═╡ 94ee248b-453f-48fc-94d6-94464a4b7e91
solve(edges, zeros(60))

# ╔═╡ 27915ffd-7fd0-44aa-9d70-5f804d7db428
ForwardDiff.gradient(x->solve(edges, x).n, zeros(60))

# ╔═╡ a29ef69c-70d0-42e5-ae81-244b4294b097
@bind seed Slider(0:10000)

# ╔═╡ b1b8f0bb-6d5e-406d-bab7-ec60c449a65a
function fullerene()
	φ = (1+√5)/2
	res = NTuple{3,Float64}[]
	for (x, y, z) in ((0.0, 1.0, 3φ), (1.0, 2 + φ, 2φ), (φ, 2.0, 2φ + 1.0))
		for (α, β, γ) in ((x,y,z), (y,z,x), (z,x,y))
			for loc in ((α,β,γ), (α,β,-γ), (α,-β,γ), (α,-β,-γ), (-α,β,γ), (-α,β,-γ), (-α,-β,γ), (-α,-β,-γ))
				if loc ∉ res
					push!(res, loc)
				end
			end
		end
	end
	return res
end;

# ╔═╡ 3c025ce2-7322-4ff1-97c2-d6aa9d3f8077
normt(x, y) = sqrt(sum(abs2, x .- y))

# ╔═╡ 6c416c44-caa6-4da3-b7e5-b790372a2b2d
let
	θ2 = 1/4*π
	ϕ2 = 1/5*π
	n = 60
	c60_xy = fullerene()
	c60_edges = [(i,j) for i=1:n, j=1:n if i<j && normt(c60_xy[i], c60_xy[j]) < 2.5]
	mask = rand(Bool, 60)
	Compose.set_default_graphic_size(12cm, 10cm)
	cam_position = SVector(0.0, 0.0, 0.5)
	rot = RotY(θ2)*RotX(ϕ2)
	cam_transform = PerspectiveMap() ∘ inv(AffineMap(rot, rot*cam_position))
	Nx = Ny = Nz = 4
	tb = textstyle(:default, fontsize(7))
	nb1 = nodestyle(:circle, fill("#FFCA22"), stroke("transparent"); r=0.01)
	nb2 = nodestyle(:circle, fill("black"), stroke("transparent"); r=0.01)
	eb = bondstyle(:default; r=0.01)
	x(i,j,k) = cam_transform(SVector(i,j,k) .* 0.03).data
	
	fig = canvas() do
		for (s, (i,j,k)) in enumerate(c60_xy)
			(mask[s] ? nb1 : nb2) >> x(i,j,k)
		end
		for (i, j) in c60_edges
			eb >> (x(c60_xy[i]...), x(c60_xy[j]...))
		end
	end
	Compose.compose(context(0.5,0.4, 1.0, 1.2), fig)
end

# ╔═╡ 60ae8afb-e005-4df1-bbcb-e601d3162bb7
function viz_config(g, config; r=0.3/sqrt(nv(g)+1), show_edgeindex=false,
        node_fontsize=100pt/sqrt(nv(g)+1),
        edge_fontsize=200pt/sqrt(nv(g)+1),
        labels=1:nv(g),
        locs=GraphPlot.spring_layout(g),
        linecolor="skyblue",
        node_edgecolor="black",
        node_facecolor1="white",
		node_facecolor2="green",
		scales=ones(nv(g)),
		shade_color="#999999",
		shade_offset = (0.01, 0.01),
    )
	xs, ys = locs
	xs .*=  scales
	ys .*=  scales
	midscale = median(scales)
	xs .= (xs .- minimum(xs)) ./ (maximum(xs)-minimum(xs))
	ys .= (ys .- minimum(ys)) ./ (maximum(ys)-minimum(ys))
    nt = nv(g)
    eb = bondstyle(:default, linewidth(5mm/sqrt(nt+1)), stroke(linecolor))
    tb1 = textstyle(:default, fontsize(node_fontsize))
    tb2 = textstyle(:default, fontsize(edge_fontsize))
	frac = -0.05
    img1 = Compose.compose(Compose.context(r, r, 1-2r, 1-2r), canvas() do
        for i=1:nt
            for j=i+1:nt
                if has_edge(g, i, j) && scales[i] > midscale
                    eb >> ((xs[i], ys[i]), (xs[j], ys[j]))
                end
            end
        end
		end
		)
	img2 = Compose.compose(Compose.context(r, r, 1-2r, 1-2r), canvas() do	
        for (x, y, s, c) in zip(xs, ys, scales, config)
			nb = if c == 1
				nodestyle(:default, fill(node_facecolor1), stroke(node_edgecolor), linewidth(2mm/sqrt(nt+1)); r=r*s)
			else
				nodestyle(:default, fill(node_facecolor2), stroke(node_edgecolor), linewidth(2mm/sqrt(nt+1)); r=r*s)
			end
			shade=nodestyle(:default, fill(shade_color), stroke(shade_color), linewidth(2mm/sqrt(nt+1)); r=r*s)
			if s > midscale
            	nb >> (x, y)
            	shade >> ((x, y) .+ shade_offset .* s)
			end
        end
    end)
	
    img3 = Compose.compose(Compose.context(r, r, 1-2r, 1-2r), canvas() do
			
        for i=1:nt
            for j=i+1:nt
                if has_edge(g, i, j) && scales[i] <= midscale
                    eb >> ((xs[i], ys[i]), (xs[j], ys[j]))
                end
            end
        end
		end
		)
	img4 = Compose.compose(Compose.context(r, r, 1-2r, 1-2r), canvas() do	
        for (x, y, s, c) in zip(xs, ys, scales, config)
			nb = if c == 1
				nodestyle(:default, fill(node_facecolor1), stroke(node_edgecolor), linewidth(1mm/sqrt(nt+1)); r=r*s)
			else
				nodestyle(:default, fill(node_facecolor2), stroke(node_edgecolor), linewidth(1mm/sqrt(nt+1)); r=r*s)
			end
			shade=nodestyle(:default, fill(shade_color), stroke(shade_color), linewidth(2mm/sqrt(nt+1)); r=r*s)
			if s <= midscale
            	nb >> (x, y)
            	shade >> ((x, y) .+ shade_offset .* s)
			end
        end
    end)
	Compose.compose(context(), 
		(context(), img2),
		(context(), img1),
		(context(), img4),
		(context(), img3),
	)

end

# ╔═╡ 4f36b93d-9e38-4ac5-b713-4efd6c788e6f
img = let
	Random.seed!(seed)
	n = 60
	locs = fullerene()
	locs = map(l->Matrix(RotX(1/3*π)*RotY(1/2*π)) * [l...], locs)
	sort!(locs, by=x->x[3])
	edges = [LightGraphs.SimpleEdge(i, j) for i in 1:n, j in 1:n if (i<j) && normt(locs[i], locs[j]) < 2.5]
	config = ForwardDiff.gradient(x->solve(map(e->(e.src=>e.dst), edges), x).n, zeros(60))
	#config = rand(0:1,60)
	scales = map(l->-l[3]*0.05+0.95, locs)
	xs = getindex.(locs, 1)
	ys = getindex.(locs, 2)
	g = SimpleGraph(edges)
	CS = ("#FF8852", "#55BE9C")
	CS = ("#FFCC52", "#55BE9C")
	CS = ("#FFCC52", "#55BEFC")
	CS = ("#00CC82", "#CC6E9C")
	CS = ("#FFECEC", "#336E8C")
	CS = ("#2C4C62", "#8C2E3C")
	LC = "#646455"
	LC = "#111111"
	img = viz_config(g, config; node_facecolor1=CS[1], linecolor=LC, node_fontsize=0, node_facecolor2=CS[2], node_edgecolor="black",
	locs=(xs, ys), scales=scales, shade_color="transparent",
		shade_offset=(0.0075, 0.0075)
	)
	Compose.set_default_graphic_size(12cm, 12cm)
	Compose.compose(context(0, 0, 0.99, 0.99), img)
end

# ╔═╡ 8ae20b2a-2527-4016-8461-ef0cb9f09443
let
	fname = "fullerene"
	img |> SVG(fname * ".svg")
	run(`rsvg-convert -f pdf -o $fname.pdf $fname.svg`)
end

# ╔═╡ d24a0ec0-44ec-4f49-86be-3e5bfdcf47fe
let
	Random.seed!(seed)
	n = 220
	config = rand(0:1, n)
	g = LightGraphs.random_regular_graph(n, 3)
	img = viz_config(g, config; node_facecolor1="#FF8852", linecolor="#646455", node_fontsize=0, node_facecolor2="#55BE9C", node_edgecolor="transparent")
	#leftright(img, updown(html"""<img src="https://user-images.githubusercontent.com/6257240/109566350-bb977f00-7ab1-11eb-953f-127d7919e3e6.png" width=270px/>""", md"Wall clock time for computing the ground state energy of Ising spin glass on the chimera graph with the ``L \times L`` unit cell (``8L^2`` spins)."))
	Compose.set_default_graphic_size(12cm, 12cm)
	img
end

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
# ╠═4dac28aa-c273-454f-981e-0087e3dbd939
# ╠═b65e6868-efdc-4b34-938a-6d5fbbad1e16
# ╠═a2edc836-42fd-4723-822b-540019980726
# ╠═94ee248b-453f-48fc-94d6-94464a4b7e91
# ╠═27915ffd-7fd0-44aa-9d70-5f804d7db428
# ╟─a29ef69c-70d0-42e5-ae81-244b4294b097
# ╠═b1b8f0bb-6d5e-406d-bab7-ec60c449a65a
# ╠═939d981f-5364-46e5-90a6-cb23d7cf8f7f
# ╠═3c025ce2-7322-4ff1-97c2-d6aa9d3f8077
# ╠═4f36b93d-9e38-4ac5-b713-4efd6c788e6f
# ╠═8ae20b2a-2527-4016-8461-ef0cb9f09443
# ╠═e4ff4e4e-6584-4a7b-839c-892647365201
# ╠═6c416c44-caa6-4da3-b7e5-b790372a2b2d
# ╠═d24a0ec0-44ec-4f49-86be-3e5bfdcf47fe
# ╠═208080a7-392d-43af-aebd-3f1d8e5af325
# ╠═5953b211-2f32-4e1c-9ebc-77ce07c5bcc5
# ╠═60ae8afb-e005-4df1-bbcb-e601d3162bb7
