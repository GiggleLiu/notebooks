### A Pluto.jl notebook ###
# v0.17.2

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

# ╔═╡ c1b529a7-478d-471e-8116-c0bb2bb52a0d
using Revise, Pkg; Pkg.activate()

# ╔═╡ a858ec12-ae05-433c-8850-84630ae10896
begin
	using PlutoUI, CoordinateTransformations, Rotations, Viznet, Compose, StaticArrays
	# left right layout
	function leftright(a, b; width=600)
		HTML("""
<style>
table.nohover tr:hover td {
   background-color: white !important;
}</style>
			
<table width=$(width)px class="nohover" style="border:none">
<tr>
	<td>$(html(a))</td>
	<td>$(html(b))</td>
</tr></table>
""")
	end
	
	# up down layout
	function updown(a, b; width=nothing)
		HTML("""<table class="nohover" style="border:none" $(width === nothing ? "" : "width=$(width)px")>
<tr>
	<td>$(html(a))</td>
</tr>
<tr>
	<td>$(html(b))</td>
</tr></table>
""")
	end
	
	function highlight(str)
		HTML("""<span style="background-color:yellow">$(str)</span>""")
	end
end;

# ╔═╡ 9c3b4290-138c-406c-a980-966dd93b11d2
using SymEngine # install with `] add SymEngine`

# ╔═╡ 7209ba2c-94d9-4102-9772-68e11c9f9fa5
using OMEinsum

# ╔═╡ a5ca245f-d67b-48a5-a70e-afdc3ee840ec
using Zygote

# ╔═╡ 62e98e87-800d-4780-9fd7-70b82f2a0c3e
using CUDA

# ╔═╡ 70cd5e7d-9200-444f-b974-0238c74c6013
using OMEinsumContractionOrders

# ╔═╡ ada57188-e70e-4ade-83d3-777192169b2b
using TropicalNumbers

# ╔═╡ bb3d0c31-787c-44d6-a415-4f85a20eac5d
html"<button onclick='present();'>present</button>"

# ╔═╡ b526ffe6-5c4a-450e-99ba-fc44fc4f6f72
html"""
<script>
document.body.style.cursor = "pointer";
</script>
"""

# ╔═╡ bac7b9bb-8d9a-418e-b8d5-ffb00796dd22
md"## Your animal friends"

# ╔═╡ ab5c94d0-a69e-476c-9548-b522c5c48fc3
rand_animal(size...) = Basic.(Symbol.('🐀' .+ rand(0:60, size...)))

# ╔═╡ e08fdfd1-a521-4434-af61-bfed3e42c61a
md"This is a vector"

# ╔═╡ dddaf055-f769-4b80-8a97-ff1e76ad99a8
v1 = rand_animal(3)

# ╔═╡ b3cba07a-84cd-4796-aa01-01e4719e8ff0
md"This is a matrix"

# ╔═╡ 1aef6d34-41f0-4717-9b13-f81f19e8825f
m1 = rand_animal(3, 3)

# ╔═╡ 399afc64-039f-43da-9994-05f4cc44d633
md"This is a tensor of rank 4"

# ╔═╡ 7a9e9466-03fd-4955-8e5d-5cba4544e7c7
t1 = rand_animal(2, 2, 2, 2)

# ╔═╡ b4088ca7-fadd-4491-a15b-9648004b238c
md"## Matrix-matrix multiplication"

# ╔═╡ 39149c34-76c1-46eb-9c1c-a0e7e2e51d17
md"This is an einsum notation, defined by `@ein_str` string literal."

# ╔═╡ ae4ef95d-0a65-4a80-8abd-1a684c60b30d
code_mm = ein"ij,jk->ik"

# ╔═╡ b3240fff-de5d-475b-9875-9440e61a858a
mm_1 = rand_animal(3, 3)

# ╔═╡ 700874d1-1dbf-4ac8-91d7-523733a176e7
mm_2 = rand_animal(3, 3)

# ╔═╡ a1b9e48d-7bb5-461c-b161-aae154fcb799
md"You can call it"

# ╔═╡ 1e9ca8b5-81bc-4946-a8c4-c08afe903d8b
ein"ij,jk->ik"(mm_1, mm_2)

# ╔═╡ 31e06d5b-2efa-41a8-987f-5e75c87a3eff
md"It is equivalent to"

# ╔═╡ 7a7dd9a1-fed7-428c-899e-90cbe3b6d699
let
	mm_out = zeros(Basic, 3, 3)
	for i=1:3
		for j=1:3
			for k=1:3
				mm_out[i,k] += mm_1[i,j] * mm_2[j, k]
			end
		end
	end
	mm_out
end

# ╔═╡ ed6fa023-43b1-4f5c-be3f-f93f16dfe285
md"It can also be constructed as"

# ╔═╡ 778d16ea-67f1-4903-9767-28682b99527c
EinCode([['i', 'j'], ['j', 'k']], ['i', 'k'])

# ╔═╡ 5cf51764-3138-4fd9-8bcd-c399ba4306f7
md"or"

# ╔═╡ f23efdd0-47dd-49b8-88c5-3debd35cd55a
EinCode([[1, 2], [2, 3]], [1, 3])

# ╔═╡ 005fddc2-1378-426c-b4b7-80a09baf8ec7
let n = Basic(:n)
	# NOTE: `flop` counts the number of iterations!
	flop(code_mm, Dict('i'=>n, 'j'=>n, 'k'=>n))
end

# ╔═╡ dfcb8259-7f6f-40b1-ae75-15eaaabe7279
md"or, for convenience"

# ╔═╡ db6ba19e-2425-4ddc-a67a-d68ebadfd0b2
flop(code_mm, uniformsize(code_mm, Basic(:n)))

# ╔═╡ 3dbb855a-4045-4433-93c8-43997161807b
md"##  This is summation"

# ╔═╡ fbb23079-6b77-412e-8648-e0bf10fabab5
sum_1 = rand_animal(3)

# ╔═╡ 7ae59add-d3d4-42bf-83fd-989c2848225a
ein"i->"(sum_1)

# ╔═╡ 14ad6f1c-cabf-47c8-ad39-81bcd92bf840
sum_2 = rand_animal(2, 2, 3)

# ╔═╡ d14e10af-5179-476d-ab8a-63d1837e866d
ein"ijk->k"(sum_2)

# ╔═╡ 61545cc3-f486-4dea-83f8-448eb23646fb
flop(ein"ijk->k", uniformsize(ein"ijk->k", Basic(:n)))

# ╔═╡ 2790a9b1-f706-4167-8b7c-8d8432b2511a
md"## Repeating a vector"

# ╔═╡ f8cc753b-8ba0-4212-80ef-1ced60866389
rp_1 = rand_animal(3)

# ╔═╡ 051de93e-f189-4a75-bd82-ebd646eb7eda
ein"i->ij"(rp_1; size_info=Dict('j'=>4))

# ╔═╡ 68e643ba-f852-4b9e-9799-13f193df3bfe
let
	rp_out = zeros(Basic, 3, 4)
	for i=1:3
		for j=1:4
			rp_out[i, j] += rp_1[i]
		end
	end
	rp_out
end

# ╔═╡ d00afc01-ecb0-454e-8cac-342d5b75101e
md"## Star contraction"

# ╔═╡ b7b777bc-3f38-4a4f-a438-454187602962
star_1 = rand_animal(2, 2)

# ╔═╡ e4ffef2b-9b19-40ed-8547-63673a3ea649
star_2 = rand_animal(2, 2)

# ╔═╡ 6c2955de-be35-4fd9-a5c8-ef430ed0013a
star_3 = rand_animal(2, 2)

# ╔═╡ 579222f1-da51-4afa-9087-2c0799bbe6ad
ein"ai, aj, ak->ijk"(star_1, star_2, star_3)

# ╔═╡ b22796e6-d163-4bbb-b4d4-bbaf57e63f8f
let
	star_out = zeros(Basic, 2, 2, 2)
	for i=1:2
		for j=1:2
			for k=1:2
				for a=1:2
					star_out[i, j, k] += star_1[a,i] * star_2[a,j] * star_3[a,k]
				end
			end
		end
	end
	star_out
end

# ╔═╡ 955d364d-5225-4833-95d5-18265e319523
md"## Automatic differentiation"

# ╔═╡ d7c25437-2bc9-473d-846e-d9f8a4ca7d9c
a, b = randn(2, 2), randn(2);

# ╔═╡ 1f57a827-c57b-4d1e-ac73-e41d6206052a
Zygote.gradient(x->ein"i->"(ein"ij,j->i"(x, b))[], a)[1]

# ╔═╡ 5dbedaab-48d8-4abe-8046-eee0e4cd4093
let
	A, B, C = randn(2,2,2,2), randn(2,2,2,2,2), randn(2,2,2,2,2)
	size_dict = uniformsize(ein"ijkl,lmkcd,asedf->dfas", 2)
	O = ein"ijkl,lmkcd,asedf->dfas"(A, B, C; size_info=size_dict)

	Ō = randn(2,2,2,2)
	# exchange input/output labels and tensors
	Ā = ein"dfas,lmkcd,asedf->ijkl"(Ō, B, C; size_info=size_dict)
	B̄ = ein"ijkl,dfas,asedf->lmkcd"(A, Ō, C; size_info=size_dict)
	C̄ = ein"ijkl,lmkcd,dfas->asedf"(A, B, Ō; size_info=size_dict)
	Ā, B̄, C̄
end;

# ╔═╡ 015b2b71-b6aa-4b6e-b6bf-75776c0f3c90
md"## Speed up your code with GPU"

# ╔═╡ 7f8a6c77-75a6-4063-9fc6-409397bb6306
md"""
* step 1: import CUDA library.
* step 2: upload your array to GPU with `CuArray` function.
"""

# ╔═╡ dcc74a07-e7e4-48a5-9942-4f7a9714db9b
let
	cuarr1, cuarr2 = CuArray(randn(2, 2)), CuArray(randn(2))
	result = ein"ij,j->i"(cuarr1, cuarr2)
	typeof(result)
end

# ╔═╡ f3b5656d-5608-4325-a884-65bceddb4963
md"""## Summary
* Einsum can be defined as: iterating over unique indices, accumulate product of corresponding input tensor elements to the output tensor.
* Einsum's representation power
    * ein"ij,jk->ik" is matrix multiplication
    * ein"i->" and ein"ijk->k" is summation
    * ein"i->ij" is repeating axis
    * ein"ai,aj,ak->ijk" is a star contraction
* The time complexity of an einsum notation is ``O(n^{(\#~ of ~ unique ~ labels)})``
* Features in OMEinsum
    * Automatic differentiation
    * GPU
    * Generic programming
"""

# ╔═╡ a9d59880-0cc1-441b-9329-f58d39c23f16
md"# Contraction order matters"

# ╔═╡ 3c50281e-ae20-40ea-85ea-dcccfbd2756e
md"## Multiplying a sequence of matrices"

# ╔═╡ 7f6c0836-7c9e-4726-9512-dc2e5ad86e1f
code_seq_1 = ein"ij,jk,kl,lm->im"

# ╔═╡ 0efbdbad-a1e1-4aab-a9e7-143a05c934eb
seq_1 = rand_animal(2,2)

# ╔═╡ e2275f15-a57b-4a74-9826-2c723aad096c
seq_2 = rand_animal(2,2)

# ╔═╡ ebfe52bb-76c1-4b39-8f8c-0dec56cc80f3
seq_3 = rand_animal(2,2)

# ╔═╡ 744576a6-6d31-48a6-a71c-cb6246b76bef
seq_4 = rand_animal(2,2)

# ╔═╡ e92b78f1-8ac6-421f-930b-e092402845cd
ein"ij,jk,kl,lm->im"(seq_1, seq_2, seq_3, seq_4)

# ╔═╡ 623551e9-48bb-4976-920e-12396a596975
flop(code_seq_1, uniformsize(code_seq_1, Basic(:n)))

# ╔═╡ d4a3557a-9de2-4500-b09e-5b3c249906fc
code_seq_2 = ein"(ij,jk),(kl,lm)->im"

# ╔═╡ 5ff2bdd9-a708-4ee3-92fa-b8fbae6bd3e3
flop(code_seq_2, uniformsize(code_seq_2, Basic(:n)))

# ╔═╡ d4040b3f-5595-4fe9-9356-6438017820fa
md"[Song Shan Lake Spring School challenge](https://github.com/QuantumBFS/SSSS)"

# ╔═╡ 484b3561-6ed3-40bf-b94d-526f508c4755
md"## The Song Shan Lake Spring School (SSSS) Challenge"

# ╔═╡ ca083761-51d9-4144-ad79-5281e3ca4522
md"In 2019, Lei Wang, Pan Zhang, Roger and me released a challenge in the Song Shan Lake Spring School, the one gives the largest number of solutions to the challenge quiz can take a macbook home ([@LinuxDaFaHao](https://github.com/LinuxDaFaHao)). Students submitted many [solutions to the problem](https://github.com/QuantumBFS/SSSS/blob/master/Challenge.md). The second part of the quiz is"

# ╔═╡ 87613668-8585-49fe-b43a-8bb4be430287
# returns atom locations
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

# ╔═╡ 8482d4d2-b377-4d36-a395-5db79ccb136c
md"""
θ = $(@bind θ2 Slider(0.0:0.01:π; default=0.5))

ϕ = $(@bind ϕ2 Slider(0.0:0.01:2π; default=0.3))
"""

# ╔═╡ dec4a621-8d5c-4c20-a9ae-d7f765852950
let
	tb = textstyle(:default)
	Compose.set_default_graphic_size(14cm, 8cm)
	cam_position = SVector(0.0, 0.0, 0.5)
	rot = RotY(θ2)*RotX(ϕ2)
	cam_transform = PerspectiveMap() ∘ inv(AffineMap(rot, rot*cam_position))
	Nx = Ny = Nz = 4
	nb = nodestyle(:circle; r=0.01)
	eb = bondstyle(:default; r=0.01)
	x(i,j,k) = cam_transform(SVector(i,j,k) .* 0.03).data
	fl = fullerene()
	fig = canvas() do
		for (i,j,k) in fl
			nb >> x(i,j,k)
			for (i2,j2,k2) in fl
				(i2-i)^2+(j2-j)^2+(k2-k)^2 < 5.0 && eb >> (x(i,j,k), x(i2,j2,k2))
			end
		end
		tb >> ((0.4, 0.2), "60 vertices\n90 edges")
		nb >> (0.4, -0.1)
		tb >> ((0.55, -0.1), "Ising spin (s=±1)")
		eb >> ((0.37, -0.05), (0.43, -0.05))
		tb >> ((0.54, -0.05), "AFM coupling")
	end
	img = Compose.compose(Compose.context(0.3,0.5, 1.2/1.4, 1.5), fig)
	img
end

# ╔═╡ 7288ba08-f845-428e-a02c-cccfe4d160cb
md"""
In the $(highlight("Buckyball")) structure shown in the figure, we attach an ising spin ``s_i=\pm 1`` on each vertex. The neighboring spins interact with an $(highlight("anti-ferromagnetic")) coupling of unit strength. Count the $(highlight("degeneracy")) of configurations that minimizes the energy
```math
E(\{s_1,s_2,\ldots,s_n\}) = \sum_{i,j \in edges}s_i s_j
```
"""

# ╔═╡ 738e9ed9-c6b3-4b90-8475-1923a3485447
md"## Step 2: Find a good contraction order"

# ╔═╡ 004ebf76-0bce-4e4d-99f6-ef60816ecc07
c60_xy = fullerene();

# ╔═╡ f3ad225a-0b42-460e-a82c-b799d5da8f13
c60_edges = [[i,j] for (i,(i2,j2,k2)) in enumerate(c60_xy), (j,(i1,j1,k1)) in enumerate(c60_xy) if i<j && (i2-i1)^2+(j2-j1)^2+(k2-k1)^2 < 5.0];

# ╔═╡ 21c55a1d-018f-480a-91e6-c46f3a516fe3
c60_code = EinCode(c60_edges, Int[])

# ╔═╡ ce134b9f-d4be-4728-b883-01b1c2ec9158
flop(c60_code, uniformsize(c60_code, Basic(:n)))

# ╔═╡ 41f07437-6f81-4da7-b4f2-037d9f36daff
c60_optcode = optimize_code(c60_code, uniformsize(c60_code, 2), TreeSA())

# ╔═╡ 647b241c-b08c-4a0c-9b40-c0398ad6b6ca
c60_elimination_order = let
	vertex_elimination_order(code::SlicedEinsum) = vertex_elimination_order(code.eins)
	vertex_elimination_order(code::NestedEinsum) = vertex_elimination_order!(code, labeltype(code)[])
	function vertex_elimination_order!(code, eliminated_vertices)
		OMEinsum.isleaf(code) && return eliminated_vertices
		for arg in code.args
			vertex_elimination_order!(arg, eliminated_vertices)
		end
		append!(eliminated_vertices, setdiff(vcat(getixsv(code.eins)...), getiyv(code.eins)))
		return eliminated_vertices
	end
	vertex_elimination_order(c60_optcode)
end

# ╔═╡ 0739e08b-0873-415e-a77b-07b83eae9fb6
md"contraction step = $(@bind nstep_c60 Slider(0:60; show_value=true, default=0))"

# ╔═╡ 739e3777-b8fc-4cd4-a92e-960f355a74ad
md"The resulting contraction order produces time complexity = $(flop(c60_optcode, uniformsize(c60_optcode, Basic(:n))))"

# ╔═╡ 78da91e7-1b76-46df-ae9a-6d611c2638e1
let
	θ2 = 0.5
	ϕ2 = 0.8
	mask = zeros(Bool,length(c60_elimination_order))
	mask[c60_elimination_order[1:nstep_c60]] .= true
	Compose.set_default_graphic_size(12cm, 12cm)
	cam_position = SVector(0.0, 0.0, 0.5)
	rot = RotY(θ2)*RotX(ϕ2)
	cam_transform = PerspectiveMap() ∘ inv(AffineMap(rot, rot*cam_position))
	Nx = Ny = Nz = 4
	tb = textstyle(:default)
	nb1 = nodestyle(:circle, fill("red"); r=0.01)
	nb2 = nodestyle(:circle, fill("white"), stroke("black"); r=0.01)
	eb = bondstyle(:default; r=0.01)
	x(i,j,k) = cam_transform(SVector(i,j,k) .* 0.03).data
	
	fig = canvas() do
		for (s, (i,j,k)) in enumerate(c60_xy)
			(mask[s] ? nb1 : nb2) >> x(i,j,k)
		end
		for (i, j) in c60_edges
			eb >> (x(c60_xy[i]...), x(c60_xy[j]...))
		end
		nb1 >> (-0.1, 0.45)
		tb >> ((-0.0, 0.45), "contracted")
		nb2 >> (-0.1, 0.50)
		tb >> ((-0.0, 0.50), "remaining")
	end
	Compose.compose(Compose.context(0.5,0.35, 1.0, 1.0), fig)
end

# ╔═╡ 174b14a9-91fe-4c6e-861d-b669fe791138
Z = c60_optcode([(J = 1.0; β = 1.0; expJ = exp(β*J); [1/expJ expJ; expJ 1/expJ]) for i=1:90]...)[]

# ╔═╡ 63001407-5b10-4b0d-b4ec-98138d24b927
log(Z) / 60

# ╔═╡ 3ce3c5ca-9402-4af8-bda1-a9c320932f0d
c60_optcode([(J = 1.0; CountingTropical.([-J J; J -J])) for i=1:90]...)

# ╔═╡ 174a1bc3-42f6-4d29-87b5-d648feb91924
md"# Summary
* Einsum is expressive, it can represent `matrix multiplication`, `batched matrix multiplication`, `trace`, `summation`, `repeating` and `tensor network`s.
* The algorithmic complexity of the einsum notation is ``n^{(\#~of~labels)}``.
* Contraction order matters.
"

# ╔═╡ Cell order:
# ╠═bb3d0c31-787c-44d6-a415-4f85a20eac5d
# ╠═b526ffe6-5c4a-450e-99ba-fc44fc4f6f72
# ╠═c1b529a7-478d-471e-8116-c0bb2bb52a0d
# ╠═a858ec12-ae05-433c-8850-84630ae10896
# ╟─bac7b9bb-8d9a-418e-b8d5-ffb00796dd22
# ╠═9c3b4290-138c-406c-a980-966dd93b11d2
# ╠═ab5c94d0-a69e-476c-9548-b522c5c48fc3
# ╟─e08fdfd1-a521-4434-af61-bfed3e42c61a
# ╠═dddaf055-f769-4b80-8a97-ff1e76ad99a8
# ╟─b3cba07a-84cd-4796-aa01-01e4719e8ff0
# ╠═1aef6d34-41f0-4717-9b13-f81f19e8825f
# ╟─399afc64-039f-43da-9994-05f4cc44d633
# ╠═7a9e9466-03fd-4955-8e5d-5cba4544e7c7
# ╟─b4088ca7-fadd-4491-a15b-9648004b238c
# ╠═7209ba2c-94d9-4102-9772-68e11c9f9fa5
# ╟─39149c34-76c1-46eb-9c1c-a0e7e2e51d17
# ╠═ae4ef95d-0a65-4a80-8abd-1a684c60b30d
# ╠═b3240fff-de5d-475b-9875-9440e61a858a
# ╠═700874d1-1dbf-4ac8-91d7-523733a176e7
# ╟─a1b9e48d-7bb5-461c-b161-aae154fcb799
# ╠═1e9ca8b5-81bc-4946-a8c4-c08afe903d8b
# ╟─31e06d5b-2efa-41a8-987f-5e75c87a3eff
# ╠═7a7dd9a1-fed7-428c-899e-90cbe3b6d699
# ╟─ed6fa023-43b1-4f5c-be3f-f93f16dfe285
# ╠═778d16ea-67f1-4903-9767-28682b99527c
# ╟─5cf51764-3138-4fd9-8bcd-c399ba4306f7
# ╠═f23efdd0-47dd-49b8-88c5-3debd35cd55a
# ╠═005fddc2-1378-426c-b4b7-80a09baf8ec7
# ╟─dfcb8259-7f6f-40b1-ae75-15eaaabe7279
# ╠═db6ba19e-2425-4ddc-a67a-d68ebadfd0b2
# ╟─3dbb855a-4045-4433-93c8-43997161807b
# ╠═fbb23079-6b77-412e-8648-e0bf10fabab5
# ╠═7ae59add-d3d4-42bf-83fd-989c2848225a
# ╠═14ad6f1c-cabf-47c8-ad39-81bcd92bf840
# ╠═d14e10af-5179-476d-ab8a-63d1837e866d
# ╠═61545cc3-f486-4dea-83f8-448eb23646fb
# ╟─2790a9b1-f706-4167-8b7c-8d8432b2511a
# ╠═f8cc753b-8ba0-4212-80ef-1ced60866389
# ╠═051de93e-f189-4a75-bd82-ebd646eb7eda
# ╠═68e643ba-f852-4b9e-9799-13f193df3bfe
# ╟─d00afc01-ecb0-454e-8cac-342d5b75101e
# ╠═b7b777bc-3f38-4a4f-a438-454187602962
# ╠═e4ffef2b-9b19-40ed-8547-63673a3ea649
# ╠═6c2955de-be35-4fd9-a5c8-ef430ed0013a
# ╠═579222f1-da51-4afa-9087-2c0799bbe6ad
# ╠═b22796e6-d163-4bbb-b4d4-bbaf57e63f8f
# ╟─955d364d-5225-4833-95d5-18265e319523
# ╠═a5ca245f-d67b-48a5-a70e-afdc3ee840ec
# ╠═d7c25437-2bc9-473d-846e-d9f8a4ca7d9c
# ╠═1f57a827-c57b-4d1e-ac73-e41d6206052a
# ╠═5dbedaab-48d8-4abe-8046-eee0e4cd4093
# ╟─015b2b71-b6aa-4b6e-b6bf-75776c0f3c90
# ╟─7f8a6c77-75a6-4063-9fc6-409397bb6306
# ╠═62e98e87-800d-4780-9fd7-70b82f2a0c3e
# ╠═dcc74a07-e7e4-48a5-9942-4f7a9714db9b
# ╟─f3b5656d-5608-4325-a884-65bceddb4963
# ╟─a9d59880-0cc1-441b-9329-f58d39c23f16
# ╟─3c50281e-ae20-40ea-85ea-dcccfbd2756e
# ╠═7f6c0836-7c9e-4726-9512-dc2e5ad86e1f
# ╠═0efbdbad-a1e1-4aab-a9e7-143a05c934eb
# ╠═e2275f15-a57b-4a74-9826-2c723aad096c
# ╠═ebfe52bb-76c1-4b39-8f8c-0dec56cc80f3
# ╠═744576a6-6d31-48a6-a71c-cb6246b76bef
# ╠═e92b78f1-8ac6-421f-930b-e092402845cd
# ╠═623551e9-48bb-4976-920e-12396a596975
# ╠═d4a3557a-9de2-4500-b09e-5b3c249906fc
# ╠═5ff2bdd9-a708-4ee3-92fa-b8fbae6bd3e3
# ╟─d4040b3f-5595-4fe9-9356-6438017820fa
# ╟─484b3561-6ed3-40bf-b94d-526f508c4755
# ╟─ca083761-51d9-4144-ad79-5281e3ca4522
# ╠═87613668-8585-49fe-b43a-8bb4be430287
# ╟─dec4a621-8d5c-4c20-a9ae-d7f765852950
# ╟─8482d4d2-b377-4d36-a395-5db79ccb136c
# ╟─7288ba08-f845-428e-a02c-cccfe4d160cb
# ╟─738e9ed9-c6b3-4b90-8475-1923a3485447
# ╠═004ebf76-0bce-4e4d-99f6-ef60816ecc07
# ╠═f3ad225a-0b42-460e-a82c-b799d5da8f13
# ╠═21c55a1d-018f-480a-91e6-c46f3a516fe3
# ╠═70cd5e7d-9200-444f-b974-0238c74c6013
# ╠═ce134b9f-d4be-4728-b883-01b1c2ec9158
# ╠═41f07437-6f81-4da7-b4f2-037d9f36daff
# ╠═647b241c-b08c-4a0c-9b40-c0398ad6b6ca
# ╟─0739e08b-0873-415e-a77b-07b83eae9fb6
# ╟─739e3777-b8fc-4cd4-a92e-960f355a74ad
# ╠═78da91e7-1b76-46df-ae9a-6d611c2638e1
# ╠═174b14a9-91fe-4c6e-861d-b669fe791138
# ╠═63001407-5b10-4b0d-b4ec-98138d24b927
# ╠═ada57188-e70e-4ade-83d3-777192169b2b
# ╠═3ce3c5ca-9402-4af8-bda1-a9c320932f0d
# ╟─174a1bc3-42f6-4d29-87b5-d648feb91924
