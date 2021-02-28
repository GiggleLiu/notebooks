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
using PlutoUI, TropicalNumbers

# ╔═╡ 5e3666e2-7961-11eb-2b2b-47737752159c
module SymTropical
	using Latexify, Markdown
	using TropicalNumbers
	using SymEngine
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
	Base.zero(::Type{Tropical{Basic}}) = Tropical(-Basic(:∞))
end

# ╔═╡ fa139cd8-7964-11eb-1c36-eb91c2f7f23c
using .SymTropical: Basic

# ╔═╡ 31f6ec44-7878-11eb-3d0e-bb365f592a0e
using Viznet, Compose

# ╔═╡ ed4600c8-788c-11eb-297c-29316dfcdec1
using CoordinateTransformations, StaticArrays, Rotations

# ╔═╡ 5143e258-7852-11eb-05fa-9df9234c5619
begin
	using LoopVectorization
	LoopVectorization.check_args(::Type{T}, ::Type{T}) where T<:Tropical = true
	LoopVectorization.check_type(::Type{Tropical{T}}) where {T} = LoopVectorization.check_type(T)
end

# ╔═╡ e45289b2-786d-11eb-3ed0-f1034789f9f2
begin
	using VectorizationBase
	using VectorizationBase: OffsetPrecalc, StaticBool, Bit, static, NativeTypes, Index, gep_quote, VectorIndex, AbstractSIMDVector, StridedPointer

	@inline function VectorizationBase.stridedpointer(A::AbstractArray{T}) where {T <: Tropical}
    stridedpointer(VectorizationBase.memory_reference(A), VectorizationBase.contiguous_axis(A),
        VectorizationBase.contiguous_batch_size(A), VectorizationBase.val_stride_rank(A),
        VectorizationBase.bytestrides(A), VectorizationBase.offsets(A))
end

	@inline function VectorizationBase.stridedpointer(
		ptr::Ptr{T}, ::StaticInt{C}, ::StaticInt{B}, ::Val{R}, strd::X, offsets::O
	) where {T<:Tropical,C,B,R,N,X<:Tuple{Vararg{Integer,N}},O<:Tuple{Vararg{Integer,N}}}
		VectorizationBase.StridedPointer{T,N,C,B,R,X,O}(ptr, strd, offsets)
	end
	
	# `vload` interfaces
	@inline function VectorizationBase.vload(ptr::Ptr{Tropical{T}}, i::I, m::Mask, a::A, si::StaticInt{RS}) where {A <: StaticBool, T <: NativeTypes, I <: Index, RS}
		Tropical(vload(Ptr{T}(ptr), i, m, a, si))
	end
	@inline function VectorizationBase.vload(ptr::Ptr{Tropical{T}}, i::I, a::A, si::StaticInt{RS}) where {A <: StaticBool, T <: NativeTypes, I <: Index, RS}
		Tropical(vload(Ptr{T}(ptr), i, a, si))
	end

	# `vstore!` interfaces
	@inline VectorizationBase.vstore!(ptr::StridedPointer{T}, v::T) where {T<:Tropical} = vstore!(ptr, content(v))
	@inline function VectorizationBase.vstore!(
		ptr::Ptr{Tropical{T}}, v::Tropical{Vec{N,T}}, i::VectorIndex{W}, m::AbstractSIMDVector{W}, a::A, s::S, nt::NT, si::StaticInt{RS}) where {T,W,S<:StaticBool,A<:StaticBool,NT<:StaticBool,RS,N}
		vstore!(convert(Ptr{T}, ptr), content(v), i, m, a, s, nt, si)
	end
	@inline function VectorizationBase.vstore!(
		ptr::Ptr{Tropical{T}}, v::Tropical{Vec{N,T}}, m::AbstractSIMDVector{W}, a::A, s::S, nt::NT, si::StaticInt{RS}) where {T,W,S<:StaticBool,A<:StaticBool,NT<:StaticBool,RS,N}
		vstore!(convert(Ptr{T}, ptr), content(v), m, a, s, nt, si)
	end
end

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

# ╔═╡ 1ba138d2-7963-11eb-1622-49c797062a8e
ising_bondtensor(Tropical{Basic}, Basic(:J))

# ╔═╡ 460adcbe-7964-11eb-261c-bdb704a573f3
ising_vertextensor(Tropical{Basic}, 2, Basic(:h))

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
log(1+0.01)

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

# ╔═╡ 2dc7a6d4-7870-11eb-2361-4d11a77da7b4
md"""
However, this motivate use to develop a reliable BLAS for tropical numbers to speed up these applications!
"""

# ╔═╡ 8d24b3aa-7853-11eb-0be4-23088fd5e70a
md"## Tropical GEMM"

# ╔═╡ bd1fb060-786b-11eb-076a-998aee8fa485
md"A short introduction (or link) to SIMD, multi-threading and tiling. The theoretical lower bound of tropical GEMM computing time is `2 x (the corresponding regular GEMM)` because the max `avx` instruction is two times slower."

# ╔═╡ 76227a6c-7870-11eb-32b3-8392e158059b
md"""
#### Why LoopVectorization and Octavian?

Because it is fast, we do not need to care about technical details of tiling et. al.
* tell `@avx` macro this type is compatible SIMD,
    * `LoopVectorization.check_args`
    * `LoopVectorization.check_type`

* vectorized data loading and storing
    * `VectorizationBase.stridedpointer`
    * `VectorizationBase.vload`
    * `VectorizationBase.vstore!`

* vectorized operations
    * `VectorizationBase._vzero`
    * `VectorizationBase.fma`
    * `VectorizationBase.vbroadcast`
    * `VectorizationBase.similar_no_offset`

* some ugly patches
    * `Base.promote`
"""

# ╔═╡ 82af0af2-786b-11eb-3a49-97519a15a851
md"#### Tell `@avx` this type is compatible SIMD"

# ╔═╡ 8167bf86-7852-11eb-0201-1996d24d3015
md"""
The `@avx` macro also checks the array arguments using `LoopVectorization.check_args` to try and determine
if they are compatible with the macro. If `check_args` returns false, a fall back loop annotated with `@inbounds`
and `@fastmath` is generated. Note that `VectorizationBase` provides functions such as `vadd` and `vmul` that will
ignore `@fastmath`, preserving IEEE semantics both within `@avx` and `@fastmath`.
`check_args` currently returns false for some wrapper types like `LinearAlgebra.UpperTriangular`, requiring you to
use their `parent`. Triangular loops aren't yet supported.

LoopVectorization will optimize an `@avx` loop if `check_args` on each on the indexed abstract arrays returns true.
It returns true for `AbstractArray{T}`s when `check_type(T) == true` and the array or its parent is a `StridedArray` or `AbstractRange`.

To provide support for a custom array type, ensure that `check_args` returns true, either through overloading it or subtyping `DenseArray`.
Additionally, define `pointer` and `stride` methods.
"""

# ╔═╡ c9effb62-786e-11eb-076c-6375642a6398
md"""
There are 4 interfaces of `vload`,
```julia
vload(p::AbstractStridedPointer, i::Tuple)
vload(p::AbstractStridedPointer, i::Tuple, m::Mask)
vload(p::AbstractStridedPointer, i::Unroll)
vload(p::AbstractStridedPointer, i::Unroll, m::Mask)

```
where `i` is the index and `m` is an optional mask to avoid loading where the mask is off.

`MM`s represent a vector of indices. Passing them as arguments basically means to take a slice. Adding a mask lets you avoid reading out of bounds.

The `Unroll` indices are for if we want to load multiple vectors at one time. 
"""

# ╔═╡ 00581528-78b2-11eb-0904-e9ab3fd376ce
md"#### vectorized operations"

# ╔═╡ 099384d8-7852-11eb-165a-6bbbcf097309
begin
	@inline function VectorizationBase.vbroadcast(a::Union{Val{W},StaticInt{W}}, s::Tropical{T}) where {W,T}
		Tropical(VectorizationBase.vbroadcast(a, content(s)))
	end

	@inline function VectorizationBase._vzero(::StaticInt{W}, ::Type{T}, ::StaticInt{RS}) where {W,T<:Tropical{FT},RS} where FT
		Tropical(VectorizationBase._vbroadcast(StaticInt{W}(), FT(-Inf), StaticInt{RS}()))
	end

	@inline function VectorizationBase.fma(x::Tropical{V}, y::Tropical{V}, z::Tropical{V}) where {V<:VectorizationBase.AbstractSIMD}
		Tropical(max(content(z), content(x) + content(y)))
	end

	@inline function VectorizationBase.similar_no_offset(sptr::OffsetPrecalc{T}, ptr::Ptr{Tropical{T}}) where {T}
		OffsetPrecalc(VectorizationBase.similar_no_offset(getfield(sptr, :ptr), ptr), getfield(sptr, :precalc))
	end

	# is `gep` a shorthand for "get element pointer"?
	@inline VectorizationBase.gep(ptr::Ptr{Tropical{T}}, i) where T = Ptr{Tropical{T}}(VectorizationBase.gep(Ptr{T}(ptr), i))
end

# ╔═╡ 4e473cf2-78b2-11eb-2fd3-035e416d1650
md"#### some ungly patches"

# ╔═╡ 490487e2-78b2-11eb-3e0c-f9461bc1739c
begin
	# TODO: FIX!!!!!!
	@inline function Base.promote(a::Int, b::Tropical{Vec{M,FT}}) where {M,FT}
		elem = a == 0 ? -Inf : 0.0
		Tropical(VectorizationBase.vbroadcast(StaticInt{M}(), FT(elem))), b
	end

	@inline function Base.promote(a::Int, b::Tropical{Vec{M,FT}}, c::Tropical{Vec{M,FT}}) where {M,FT}
		elem = a == 0 ? -Inf : 0.0
		Tropical(VectorizationBase.vbroadcast(StaticInt{M}(), FT(elem))), b, c
	end
end

# ╔═╡ 6856bd10-7859-11eb-0d4f-4f3342310ddb

function distance(a::AbstractArray{<:Tropical}, b::AbstractArray{<:Tropical})
    sum(abs.(content.(a) .- content.(b)))
end

@testset "matmul" begin
    for n in [4, 40]
        a = Tropical.(randn(n, n))
        b = Tropical.(randn(n, n))
        @test distance(Octavian.matmul_serial(a, b), a*b) ≈ 0
        @test distance(Octavian.matmul_serial(a, a), a*a) ≈ 0
        @test distance(Octavian.matmul(a, b), a*b) ≈ 0
    end
end

function naivemm!(o::Matrix, a::Matrix, b::Matrix)
    @assert size(a, 2) == size(b, 1) && size(o) == (size(a, 1), size(b, 2))
    for j=1:size(b, 2)
        for k=1:size(a, 2)
            for i=1:size(a, 1)
                @inbounds o[i,j] += a[i,k] * b[k,j]
            end
        end
    end
    return o
end

# ╔═╡ 94bf6aba-7871-11eb-1431-697b41de04a3
md"## Benchmarks"

# ╔═╡ 93abda52-786d-11eb-2b2e-0787a202c609
md"## Limitations"

# ╔═╡ b37ec12a-786d-11eb-286c-01ba8c3546c8
md"""With the disclaimer that those types would have to map 1-1 to native numbers for it to work well, because LoopVectorization assumes that is the case in a way critical to it's ability to optimize code. So this works for Tropical numbers, but it wouldn't (for example) Complex or ForwardDiff.Dual numbers, quarternions, or RGB colors.
I'll probably get around to making things like these work eventually using the AbstractInterpreter interface, but the "todo" list before I get there is still quite long."""

# ╔═╡ 695e405c-786d-11eb-0a6e-bb776d9626ad
md"
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

# ╔═╡ Cell order:
# ╠═c456b902-7959-11eb-03ba-dd14a2cd5758
# ╠═0faacc8a-7965-11eb-151c-2909d9c2f00e
# ╠═5e3666e2-7961-11eb-2b2b-47737752159c
# ╠═fa139cd8-7964-11eb-1c36-eb91c2f7f23c
# ╠═1ba138d2-7963-11eb-1622-49c797062a8e
# ╠═460adcbe-7964-11eb-261c-bdb704a573f3
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
# ╠═faa0b5bc-794a-11eb-186e-e5cc9e6f4b15
# ╟─2dc7a6d4-7870-11eb-2361-4d11a77da7b4
# ╟─8d24b3aa-7853-11eb-0be4-23088fd5e70a
# ╟─bd1fb060-786b-11eb-076a-998aee8fa485
# ╠═76227a6c-7870-11eb-32b3-8392e158059b
# ╟─82af0af2-786b-11eb-3a49-97519a15a851
# ╟─8167bf86-7852-11eb-0201-1996d24d3015
# ╠═5143e258-7852-11eb-05fa-9df9234c5619
# ╟─c9effb62-786e-11eb-076c-6375642a6398
# ╠═e45289b2-786d-11eb-3ed0-f1034789f9f2
# ╟─00581528-78b2-11eb-0904-e9ab3fd376ce
# ╠═099384d8-7852-11eb-165a-6bbbcf097309
# ╟─4e473cf2-78b2-11eb-2fd3-035e416d1650
# ╠═490487e2-78b2-11eb-3e0c-f9461bc1739c
# ╠═6856bd10-7859-11eb-0d4f-4f3342310ddb
# ╟─94bf6aba-7871-11eb-1431-697b41de04a3
# ╟─93abda52-786d-11eb-2b2e-0787a202c609
# ╟─b37ec12a-786d-11eb-286c-01ba8c3546c8
# ╟─695e405c-786d-11eb-0a6e-bb776d9626ad
