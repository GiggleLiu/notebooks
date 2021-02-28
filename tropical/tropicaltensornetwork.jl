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

# ╔═╡ c456b902-7959-11eb-03ba-dd14a2cd5758
using PlutoUI, TropicalNumbers, LightGraphs

# ╔═╡ 5e3666e2-7961-11eb-2b2b-47737752159c
module SymTropical
	using Latexify, Markdown
	using TropicalNumbers
	using SymEngine
	using HierarchicalBipartition
	using SimpleTensorNetworks

	function Latexify.latexraw(x::Tropical{Basic})
		return string(TropicalNumbers.content(x))
	end
	function Base.show(io::IO, mime::MIME"text/html", arr::AbstractArray{<:Tropical{Basic}})
		Base.show(io, mime, Markdown.parse("""
```math
$(latexify(arr))
```
"""))
	end
	function build_tree(N::Int, order)
		ids = collect(Any, 1:N)
		for i=1:length(order)
			ta, tb = order[i].src, order[i].dst
			ids[ta] = ContractionTree(ids[ta], ids[tb])
			ids[tb] = nothing
		end
		filter(x->x!==(nothing), ids)
	end
	Base.zero(::Type{Tropical{Basic}}) = Tropical(-Basic(:∞))
	function HierarchicalBipartition.order_greedy(tn::TensorNetwork)
		graph = SimpleTensorNetworks.tn2graph(tn)
		log2shapes = [[log2.(size(tn.tensors[i]))...] for i=1:length(tn)]
		layout = TensorNetworkLayout(graph, log2shapes)
		tc, sc, order = HierarchicalBipartition.order_greedy(layout)
		tc, sc, build_tree(length(tn), order)
	end
end

# ╔═╡ fa139cd8-7964-11eb-1c36-eb91c2f7f23c
using .SymTropical: Basic

# ╔═╡ 31f6ec44-7878-11eb-3d0e-bb365f592a0e
using Viznet, Compose

# ╔═╡ ed4600c8-788c-11eb-297c-29316dfcdec1
using CoordinateTransformations, StaticArrays, Rotations

# ╔═╡ bbdcbf6c-7984-11eb-2219-159d553fc8b5
using SimpleTensorNetworks

# ╔═╡ b8ae5d48-798c-11eb-343a-03c5fa330ea5
using Random

# ╔═╡ 0faacc8a-7965-11eb-151c-2909d9c2f00e
begin
	
	function ising_bondtensor(::Type{T}, J) where T
		e = T(J)
		e_ = T(-J)
		[e e_; e_ e]
	end
	function ising_vertextensor(::Type{T}, n::Int, h) where T
		res = zeros(T, fill(2, n)...)
		res[1] = T(h)
		res[end] = T(-h)
		return res
	end
	
	function twosat_bondtensor(::Type{T}, src::Bool, dst::Bool) where T
		res = [T(1) T(1); T(1) T(1)]
		res[Int(src)+1, Int(dst)+1] = T(-1)
		return res
	end
	
	function twosat_vertextensor(::Type{T}, n::Int) where T
		res = zeros(T, fill(2, n)...)
		res[1] = one(T)
		res[end] = one(T)
		return res
	end
	
	function potts_bondtensor(::Type{T}, ::Val{q}, J) where {T, q}
		angles = cos.(2π .* ((1:q) ./ q))
		res = zeros(T, q, q)
		for i=1:q
			for j=1:q
				res[i,j] = T(J*angles[mod1(abs(j-i), q)])
			end
		end
		res
	end

	function potts_vertextensor(::Type{T}, ::Val{q}, n::Int) where {T, q}
		res = zeros(T, fill(q, n)...)
		for i=1:q
			res[fill(i, n)...] = one(T)
		end
		res
	end


	#for F in [:ising_bondtensor, :]
	#	@eval $F(t) = 
	#end
end

# ╔═╡ be76e52a-7852-11eb-179b-afbc6efcab55
md"## Tropical algebra"

# ╔═╡ d0b54b76-7852-11eb-2398-0911380fa090
md"""
Tropical algebra is defined by replacing the usual sum and product operators for ordinary real numbers with the max and sum operators respectively~\cite{maclagan2015introduction} 

```math
\begin{align}
&a ⊕ b = \max(a, b)\\
&a ⊙ b = a + b
\end{align}
```
"""

# ╔═╡ af13e090-7852-11eb-21ae-8b94f25f1a4f
Tropical(3.0) * Tropical(2.0)

# ╔═╡ d770f232-7864-11eb-0e9a-81528e359d39
Tropical(3.0) + Tropical(2.0)

# ╔═╡ 3372871c-785b-11eb-3092-4bbc419cb788
md"One sees that $-\infty$ acts as zero element for the tropical number since  $-\infty \oplus x = x  $ and $-\infty \odot x = -\infty$. "

# ╔═╡ 2173a6cc-785b-11eb-1ab6-7fb875224dd9
zero(Tropical{Float64})

# ╔═╡ 518b7d4e-785b-11eb-3b7c-1389065b9cbd
md"
On the other hand, $0$ acts as the multiplicative identity since $0 \odot x = x$."

# ╔═╡ 2868b292-785b-11eb-015e-6b5613bd9e39
one(Tropical{Float64})

# ╔═╡ 5d16a472-785b-11eb-1b94-dd6d8f860c24
md"""
The $\oplus$ and $\odot$ operators still have commutative, associative, and distributive properties. However, since there is no additive inverse, the $\oplus$ and $\odot$ and operations define a semiring over ``\mathbb R \cup \{-\infty\}``. 
"""

# ╔═╡ 98ae0960-797d-11eb-3646-c5b7e05d3f7c
md"""The tropical $\delta$ tensor of rank $n$ and dimension $q$ is defined as
```math
δ_{s_i s_j\ldots s_n}^{n,q} = \begin{cases}
 0, & s_i = s_j =\ldots s_n\\
 -\infty, &otherwise
\end{cases}
```
where $s_i,s_j,\ldots s_n \in \{1,2,\ldots q\}$.
"""

# ╔═╡ 695e405c-786d-11eb-0a6e-bb776d9626ad
md"
## Using Tropical numbers for counting
Moreover, one can also employ the present approach to count the number of ground states at the same computational complexity of computing the ground state energy. To implement this, we further generalize the tensor element to be a tuple $(x, n)$ composed by a tropical number $x$ and an ordinary number $n$. The tropical number records the negative energy, while the ordinary number counts the number of minimal energy configurations. For tensor network contraction, we need the multiplication and addition of the tuple: $(x_1, n_1) \odot (x_2,n_2) = (x_1 + x_2, n_1\cdot n_2)$ and 
```math
\begin{equation}
    (x_1, n_1)\oplus (x_2, n_2) = \begin{cases}
 (x_1\oplus x_2, \, n_1 + n_2 ) & \text{if $x_1 = x_2$} \\
 (x_1\oplus x_2,\, n_1 ) & \text{if $x_1>x_2$} \\
 (x_1\oplus x_2,\, n_2 )& \text{if $x_1 < x_2$}
 \end{cases}.
\end{equation}
```
"

# ╔═╡ 442bcb3c-7940-11eb-18e5-d3158b74b1dc
html"""
<h2>Mapping hard problems to Tropical Tensor networks</h2>
<table style="border:none">
<tr>
	<td rowspan=4>
	<img src="https://images-na.ssl-images-amazon.com/images/I/51QttTd6JLL._SX351_BO1,204,203,200_.jpg" width=200px/>
	</td>
	<td rowspan=1 align="center">
	<big>The Nature of Computation</big><br><br>
	By <strong>Cristopher Moore</strong>
	</td>
</tr>
<tr>
	<td align="center">
	<strong>Section 5</strong>
	<br><br>Who is the hardest one of All?
	<br>NP-Completeness
	</td>
</tr>
<tr>
	<td align="center">
	<strong>Section 13</strong>
	<br><br>Counting, sampling and statistical physics
	</td>
</tr>
</table>
"""

# ╔═╡ f7208b6e-793c-11eb-0dfa-0d63752ba53e
md"""## Ising Spin glass
* hardness: NP-Complete
* counting: #P

```math
-E = \sum\limits_{i,j\in E} J_{ij} s_i s_j + \sum\limits_{i\in V} h_{i} s_i
```
```math
T_{e} = \mathcal{T}\begin{bmatrix}J & -J \\-J & J\end{bmatrix}
```

```math
T_{v}^{n} = \begin{cases}
 h, & s_i = s_j =\ldots s_n = 0\\
 -h, & s_i = s_j =\ldots s_n = 1\\
 -\infty, &otherwise
\end{cases}
```
where $s_i,s_j,\ldots s_n \in \{1,2\}$.
"""

# ╔═╡ ff77ceea-785a-11eb-2e71-6f4bc8c10881
md"""
* the exact ground state energy of Ising spin glasses on square lattice up to $32^2$ spins, on cubic lattice up to $6^3$ spins

* We obtain exact ground state energy of $\pm J$ Ising spin glass on the chimera graph of D-Wave quantum annealer of $512$ qubits in less than $100$ seconds and investigate the exact value of the residual entropy of $\pm J$ spin glasses on the chimera graph;


* The spin glass on the random graphs: our method can compute optimal solutions and count the number of solutions for spin glasses and combinatorial optimization problems on on $3$ regular random graphs up to $220$ spins, on a single GPU. This is inaccessible by existing methods.
"""

# ╔═╡ 06bbead0-793f-11eb-0dec-c549b461b9cf
md"""
## Max 2-satisfiability problem
* hardness: Polynomial
* counting: #P

A 2-satisfiability problem may be described using a Boolean expression with a special restricted form. It is a conjunction (a Boolean and operation) of clauses, where each clause is a disjunction (a Boolean or operation) of two variables or negated variables. The variables or their negations appearing in this formula are known as literals.

```math
\begin{align}
& (x_{0}\lor x_{2})\land (x_{0}\lor \lnot x_{3})\land (x_{1}\lor \lnot x_{3})\land (x_{1}\lor \lnot x_{4})\land \\
& (x_{2}\lor \lnot x_{4})\land {}(x_{0}\lor \lnot x_{5})\land (x_{1}\lor \lnot x_{5})\land (x_{2}\lor \lnot x_{5})\land \\
&(x_{3}\lor x_{6})\land (x_{4}\lor x_{6})\land (x_{5}\lor x_{6}).
\end{align}
```

The spin glass and Max $2$-SAT problem on the random graphs: our method can compute optimal solutions and count the number of solutions for spin glasses and combinatorial optimization problems on on $3$ regular random graphs up to $220$ spins, on a single GPU. This is inaccessible by existing methods.
"""

# ╔═╡ ef2d2446-793f-11eb-223a-c5fe0ed5e367
md"""
```math
E = \sum\limits_{i =1}^{|C|} C_i(\mathbf{s})
```

```math
(T_{e})_{s_i s_j} = \begin{cases}
	1,  & clause(e,s_i,s_j) ~\text{is true}\\
	-1, &otherwise
\end{cases}
```

```math
T_{v}^n = \delta^{n,q=2}
```
"""

# ╔═╡ 5f2243c4-793d-11eb-1add-392387bb559f
md"""
## Potts model
The Potts model consists of spins that are placed on a lattice; the lattice is usually taken to be a two-dimensional rectangular Euclidean lattice, but is often generalized to other dimensions or other lattices. Domb originally suggested that the spin takes one of q possible values, distributed uniformly about the circle, at angles

``\theta _{n}={\frac  {2\pi n}{q}}``,
where ``n = 0, 1, ..., q-1`` and that the interaction Hamiltonian be given by

```math
H=J\sum _{{i,j \in E}}\cos \left(\theta _{{s_{i}}}-\theta _{{s_{j}}}\right)
```
For $q=3$,
```math
T_e = J\left(\begin{matrix}1 & -1/2 & -1/2 \\ -1/2 & 1 & -1/2 \\ -1/2 & -1/2 & 1\end{matrix}\right)
```
and 

$$T_v^n=\delta^{n,q=3}$$



* we investigate ground-state energy and entropy of $3$-state Potts glasses on square lattices up to size $18\times 18$.
"""

# ╔═╡ b3d665e2-7967-11eb-13ce-f9ab0c0fd4b7
potts_vertextensor(Tropical{Basic}, Val(3), 2)

# ╔═╡ 344042b4-793d-11eb-3d6f-43eb2a4db9f4
md"""
#### Maximum independent set
```math
T_{b} = \mathcal{T}\begin{bmatrix}0 & 0 \\0 & -\infty\end{bmatrix}
```

```math
T_{v}^{n} = \begin{cases}
 0, & s_i = s_j =\ldots s_n = 0\\
 1, & s_i = s_j =\ldots s_n = 1\\
 -\infty, &otherwise
\end{cases}
```
where $s_i,s_j,\ldots s_n \in \{1,2\}$.
"""

# ╔═╡ 326f2b30-787f-11eb-0a63-6b76097d84b6
struct Cubic{T}
	grid::NTuple{3,T}
end

# ╔═╡ b1a751ea-7882-11eb-082a-4d18382cedcc
function Base.getindex(c, i, j, k)
	(i,j,k) .* c.grid
end

# ╔═╡ 9f082bce-788d-11eb-1555-fb602a81dfee
md"``R_x`` = $(@bind θ Slider(0:0.01:2π; default=0.5))"

# ╔═╡ 87531a00-7891-11eb-0072-5b9aebb4625a
md"``R_y`` = $(@bind ϕ Slider(0:0.01:2π; default=2.8))"

# ╔═╡ a5e2efee-7883-11eb-16d7-9d2769c08435
let
	Compose.set_default_graphic_size(12cm, 12cm)
	cam_position = SVector(0.0, 0.0, 0.5)
	rot = RotY(θ)*RotX(ϕ)
	cam_transform = PerspectiveMap() ∘ inv(AffineMap(rot, rot*cam_position))
	Nx = Ny = Nz = 4
	nb = nodestyle(:circle; r=0.01)
	eb = bondstyle(:default; r=0.01)
	c = Cubic((0.05, 0.05, 0.05))
	x(i,j,k) = cam_transform(SVector(c[i-Nx/2-0.5,j-Ny/2-0.5,k-Nz/2-0.5])).data
	fig = canvas() do
		for i=1:Nx, j=1:Ny, k=1:Nz
			nb >> x(i,j,k)
			i!=Nx && eb >> (x(i,j,k), x(i+1,j,k))
			j!=Ny && eb >> (x(i,j,k), x(i,j+1,k))
			k!=Nz && eb >> (x(i,j,k), x(i,j,k+1))
		end
	end
	Compose.compose(context(0.5,0.5, 1.0, 1.0), fig)
end

# ╔═╡ 5a5d4de6-7895-11eb-15c6-bda7a4342002
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
end

# ╔═╡ 1dbb9e90-78b0-11eb-2014-6dc6cfb35387
md"``R_x`` = $(@bind θ2 Slider(0:0.01:2π; default=0.5))"

# ╔═╡ 1dbc9afc-78b0-11eb-0940-2dcadf5408bb
md"``R_y`` = $(@bind ϕ2 Slider(0:0.01:2π; default=2.8))"

# ╔═╡ 9b1dc21a-7896-11eb-21f6-bfe9b4dc9ccf
let
	Compose.set_default_graphic_size(12cm, 12cm)
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
	end
	Compose.compose(context(0.5,0.5, 1.0, 1.0), fig)
end

# ╔═╡ 04e8a7da-7952-11eb-0470-d50d972083eb
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
]

# ╔═╡ 8315e3de-7984-11eb-2273-198198ee6eeb


# ╔═╡ c26b5bb6-7984-11eb-18fe-2b6a524f5c85
tn = TensorNetwork(vcat(
	[LabeledTensor(ising_vertextensor(TropicalF64, 3, 0.0), [(i,j==e.first) for (i,e) in enumerate(edges) if j ∈ e]) for j=1:60],  # vertex
	[LabeledTensor(ising_bondtensor(TropicalF64, -1.0), [(i,true),(i,false)]) for i = 1:length(edges)]
));

# ╔═╡ ae92d828-7984-11eb-31c8-8b3f9a071c24
tc, sc, order = (Random.seed!(3); order_greedy(tn))

# ╔═╡ 2b899624-798c-11eb-20c4-fd5523f7abff
length(order)

# ╔═╡ 896c80d8-798b-11eb-3319-456218053b8d
maximum(sc)

# ╔═╡ d2161642-798a-11eb-2dec-cfe6cda6af5c
SimpleTensorNetworks.contract(tn, order)

# ╔═╡ faa0b5bc-794a-11eb-186e-e5cc9e6f4b15
md"## Partition function
```math
Z = \sum\limits_{\boldsymbol{\sigma}} \prod_{i,j \in E} e^{-\beta J_{ij} \sigma_i\sigma_j}\prod_{i\in V}e^{-\beta h_i\sigma_i}
```
```math
\begin{align}
E^* &= \lim_{\beta \rightarrow \infty}-\frac{1}{\beta}\ln Z \\
&= \max\limits_{\boldsymbol{\sigma}} \sum_{i,j \in E} J_{ij} \sigma_i\sigma_j + \sum_{i\in V}h_i\sigma_i
\end{align}
```

```math
\begin{align}
\lim_{\beta\rightarrow \infty}\frac{1}{\beta} \ln (e^{\beta x} + e^{\beta y})&= x \oplus y \\
\frac{1}{\beta}\ln ( e^{\beta x} \cdot e^{\beta y}) &= x \odot y
\end{align}
```
"

# ╔═╡ 93abda52-786d-11eb-2b2e-0787a202c609
md"## Limitations"

# ╔═╡ Cell order:
# ╠═c456b902-7959-11eb-03ba-dd14a2cd5758
# ╟─0faacc8a-7965-11eb-151c-2909d9c2f00e
# ╠═5e3666e2-7961-11eb-2b2b-47737752159c
# ╠═fa139cd8-7964-11eb-1c36-eb91c2f7f23c
# ╟─be76e52a-7852-11eb-179b-afbc6efcab55
# ╟─d0b54b76-7852-11eb-2398-0911380fa090
# ╠═af13e090-7852-11eb-21ae-8b94f25f1a4f
# ╠═d770f232-7864-11eb-0e9a-81528e359d39
# ╟─3372871c-785b-11eb-3092-4bbc419cb788
# ╠═2173a6cc-785b-11eb-1ab6-7fb875224dd9
# ╟─518b7d4e-785b-11eb-3b7c-1389065b9cbd
# ╠═2868b292-785b-11eb-015e-6b5613bd9e39
# ╟─5d16a472-785b-11eb-1b94-dd6d8f860c24
# ╟─98ae0960-797d-11eb-3646-c5b7e05d3f7c
# ╟─695e405c-786d-11eb-0a6e-bb776d9626ad
# ╟─442bcb3c-7940-11eb-18e5-d3158b74b1dc
# ╟─f7208b6e-793c-11eb-0dfa-0d63752ba53e
# ╟─ff77ceea-785a-11eb-2e71-6f4bc8c10881
# ╟─06bbead0-793f-11eb-0dec-c549b461b9cf
# ╟─ef2d2446-793f-11eb-223a-c5fe0ed5e367
# ╟─5f2243c4-793d-11eb-1add-392387bb559f
# ╠═b3d665e2-7967-11eb-13ce-f9ab0c0fd4b7
# ╟─344042b4-793d-11eb-3d6f-43eb2a4db9f4
# ╠═31f6ec44-7878-11eb-3d0e-bb365f592a0e
# ╠═326f2b30-787f-11eb-0a63-6b76097d84b6
# ╠═b1a751ea-7882-11eb-082a-4d18382cedcc
# ╠═ed4600c8-788c-11eb-297c-29316dfcdec1
# ╟─9f082bce-788d-11eb-1555-fb602a81dfee
# ╟─87531a00-7891-11eb-0072-5b9aebb4625a
# ╟─a5e2efee-7883-11eb-16d7-9d2769c08435
# ╟─5a5d4de6-7895-11eb-15c6-bda7a4342002
# ╟─1dbb9e90-78b0-11eb-2014-6dc6cfb35387
# ╟─1dbc9afc-78b0-11eb-0940-2dcadf5408bb
# ╟─9b1dc21a-7896-11eb-21f6-bfe9b4dc9ccf
# ╠═04e8a7da-7952-11eb-0470-d50d972083eb
# ╠═8315e3de-7984-11eb-2273-198198ee6eeb
# ╠═bbdcbf6c-7984-11eb-2219-159d553fc8b5
# ╠═c26b5bb6-7984-11eb-18fe-2b6a524f5c85
# ╠═b8ae5d48-798c-11eb-343a-03c5fa330ea5
# ╠═ae92d828-7984-11eb-31c8-8b3f9a071c24
# ╠═2b899624-798c-11eb-20c4-fd5523f7abff
# ╠═896c80d8-798b-11eb-3319-456218053b8d
# ╠═d2161642-798a-11eb-2dec-cfe6cda6af5c
# ╟─faa0b5bc-794a-11eb-186e-e5cc9e6f4b15
# ╟─93abda52-786d-11eb-2b2e-0787a202c609
