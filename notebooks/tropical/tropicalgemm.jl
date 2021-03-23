### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

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

# ╔═╡ f0b07c72-8999-11eb-316a-293e6cee3c88
using Plots

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
end

# ╔═╡ 1ba138d2-7963-11eb-1622-49c797062a8e
ising_bondtensor(Tropical{Basic}, Basic(:J))

# ╔═╡ 460adcbe-7964-11eb-261c-bdb704a573f3
ising_vertextensor(Tropical{Basic}, 2, Basic(:h))

# ╔═╡ 8d24b3aa-7853-11eb-0be4-23088fd5e70a
md"# Tropical GEMM"

# ╔═╡ 56082ee0-898f-11eb-13fc-ab4eb456e479
md"This blog is about how to make a GEMM extension for Tropical numbers, with a close to theoretical optimal performance. It is based on
* `LoopVectorization`, for vectorizing loops (i.e. utilizing SIMD),
* and `Octavian`, a native Julia GEMM library with similar to MKL performance.
"

# ╔═╡ 31501e08-899d-11eb-0f3c-d95f668c0990
md"To appreciate the tropical GEMM better, we highly recommend readers to read this pluto notebook about [Tropical tensor networks](https://giggleliu.github.io/notebooks/tropical/tropicaltensornetwork.html)."

# ╔═╡ def82aee-898e-11eb-3b8d-2325f3709f73
md"""## Warnings before reading

The method introduced to make a BLAS extension is not garanteed to work for other user defined types. The types would have to map 1-1 to native numbers for it to work well, because LoopVectorization assumes that is the case in a way critical to it's ability to optimize code. So this works for `Tropical` numbers, but it wouldn't (for example) `Complex` or `ForwardDiff.Dual` numbers, `quarternions`, or `RGB` colors. (Chris Elrod: I'll probably get around to making things like these work eventually using the AbstractInterpreter interface, but the "todo" list before I get there is still quite long.)
"""

# ╔═╡ dba0b4f6-8993-11eb-1822-c993a037dc6b
md"## What is the goal?"

# ╔═╡ 9bbaefa0-8993-11eb-37a9-854dea2f12dd
md"""
The CPU we used for testing is `Intel(R) Core(TM) i5-10400 CPU @ 2.90GHz`. We want to sqeeze every drop of its computing power. Its theoretical serial computing power for computing a Float64 matrix multiplier is

	Serial CPU power = 2.9 GHz (CPU clock speed, we use the maximum Turbo frequency)
				  * 2 (multiplication and add can happen at the same CPU clock)
				  * 2 (number of instructions per cycle) Q: what is this?
			      * 4 (avx instruction set has a 256 with register, it can
                       crunch 4 vectorized double precision floating point
					   operations at one CPU cycle)
				= 68.8 GFLOPS

The easist way to determine a CPU's computing power is by running a matrix multiplication of size 1000 x 1000,

```julia
julia> using LinearAlgebra

julia> BLAS.vendor()  # super important to use MKL for Intel CPUs
:mkl

julia> BLAS.set_num_threads(1)

julia> @btime Octavian.matmul_serial!($(zero(A)), $A, $A);
  37.352 ms (0 allocations: 0 bytes)

julia> GFLOPS_OCT = 2 / 37.352e-3  # there are 2^9 floating point operations in total
53.54465624330691

julia> @btime LinearAlgebra.mul!($(zero(A)), $A, $A);
  33.627 ms (0 allocations: 0 bytes)

julia> GFLOPS_MKL = 2 / 33.627e-3  # there are 2^9 floating point operations in total
59.47601629642847
```
"""

# ╔═╡ bd1fb060-786b-11eb-076a-998aee8fa485
md"The theoretical computing power for tropical matrix multplication is smaller than regular floating point matrix multplication by a factor of 2, because it does not have `fma` like shortcut to do `*` and `+` in a same CPU cycle. So the best we can expect is

```julia
julia> 2/34.4
0.05813953488372093   # 58.1ms
```

Let's see how much we can approach this limit with `Octavian`.
"

# ╔═╡ 76227a6c-7870-11eb-32b3-8392e158059b
md"""
## Implementations
### Why LoopVectorization and Octavian?

It is fast and extensible. It has devided the problems into small pieces, and handled technical details of tiling et al. What people need to do is determining how to load and store data, how to do elementary arithmetic operations.
"""

# ╔═╡ ca278bd6-89a3-11eb-2388-1d50ae560b7c
md"""### Concepts
##### 1. The number types you need to handle
* `Tropical{VectorizationBase.NativeTypes}`
* `Tropical{<:VectorizationBase.Vec}`, a vector of Tropical numbers
* `Tropical{<:VectorizationBase.VecUnroll}`, a bundle of `Vec`s

Here, `NativeTypes` incldues common floating point numbers, integers, and bit types. Here we use `Tropical{<:Vec}` to present a vector (that can fit into an SIMD register) of Tropical numbers rather than something like `Vec{N, <:Tropical}` because `Vec` is finally processed by SIMD instructions, it can only contain certain `NativeTypes`.
"""

# ╔═╡ c1cc1e72-89a6-11eb-3c1d-8ba9aa0b5bb4
vec = Vec(1.0, 2.0, 3.0, 4.0)

# ╔═╡ cbb85948-89a6-11eb-092e-f1194a9774d9
vec_unroll = VecUnroll((vec, vec))

# ╔═╡ f0141896-89a6-11eb-05fe-9d140d242105
md"The reason why we need `VecUnroll` is because it is often faster to unroll a small bundle of vectorized instructions in a loop."

# ╔═╡ 64af7e68-89b1-11eb-3d43-c1a80cb69dcd
md"""
##### 2. Masks
A mask is mainly used to avoid loading/storing elements out of bounds (Q: it is correct to say out of bounds?).
When overload an interface, we often implement both the masked and the non-masked versions.
"""

# ╔═╡ 6cf7c3e6-89b1-11eb-3374-83471f184496
m = VectorizationBase.Mask{4}(0xe)

# ╔═╡ 8d9f1b94-89b1-11eb-12c2-991d640e839b
VectorizationBase.vload(pointer([1,2,3,4,5]), 2)

# ╔═╡ e84b96f2-89a3-11eb-3f92-4d4b2849b2e0
md"""
##### The interfaces that you need to implement

* vectorized data loading and storing
    * `VectorizationBase.stridedpointer`
    * `VectorizationBase.gep`, for get element pointer
    * `VectorizationBase.VectorizationBase.__vload`
    * `VectorizationBase._vstore!` and `VectorizationBase.__vstore!`

* vectorized operations
    * `VectorizationBase._vzero`
    * `VectorizationBase.zero_vecunroll`
    * `Base.fma`
    * `VectorizationBase._vbroadcast`
    * `VectorizationBase.vsum`
    * `VectorizationBase.ifelse`
    * `VectorizationBase.reduce_add`
    * `VectorizationBase.contract_add`
    * `VectorizationBase.similar_no_offset`
    * `Base.FastMath.add_fast`
    * `Base.FastMath.mul_fast`

* tell `@avx` macro this type is compatible SIMD,
    * `LoopVectorization.check_args`
    * `LoopVectorization.check_type`

* other interfaces
    * `Base.promote_rule`
    * `VectorizationBase.vecmemaybe`  Q: what is this for

Since these functions are very restrictve on types, it is easy to figure out the interfaces one needs to overwrite by try-and-error.
"""

# ╔═╡ 82af0af2-786b-11eb-3a49-97519a15a851
md"#### Step 1: Tell `@avx` Tropical numbers are is compatible with SIMD"

# ╔═╡ 8167bf86-7852-11eb-0201-1996d24d3015
md"""
The first thing is telling `@avx` macro that the `Tropical` type can utilize SIMD to avoid running into the fallback implementations.
`@avx` is a macro in `LoopVectorization` that designed to vectorize a loop automatically, it is the corner stone of `Octavian`.

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

# ╔═╡ e39c24c8-8999-11eb-0def-7be2b449dd4f
md"""
```julia
julia> @benchmark Octavian.matmul_serial!($(zero(a)), $a, $a)
BenchmarkTools.Trial: 
  memory estimate:  0 bytes
  allocs estimate:  0
  --------------
  minimum time:     63.465 ms (0.00% GC)
  median time:      63.739 ms (0.00% GC)
  mean time:        63.900 ms (0.00% GC)
  maximum time:     65.141 ms (0.00% GC)
  --------------
  samples:          79
  evals/sample:     1
```
"""

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
# ╟─0faacc8a-7965-11eb-151c-2909d9c2f00e
# ╟─5e3666e2-7961-11eb-2b2b-47737752159c
# ╠═fa139cd8-7964-11eb-1c36-eb91c2f7f23c
# ╠═1ba138d2-7963-11eb-1622-49c797062a8e
# ╠═460adcbe-7964-11eb-261c-bdb704a573f3
# ╟─8d24b3aa-7853-11eb-0be4-23088fd5e70a
# ╟─56082ee0-898f-11eb-13fc-ab4eb456e479
# ╟─31501e08-899d-11eb-0f3c-d95f668c0990
# ╟─def82aee-898e-11eb-3b8d-2325f3709f73
# ╟─dba0b4f6-8993-11eb-1822-c993a037dc6b
# ╟─9bbaefa0-8993-11eb-37a9-854dea2f12dd
# ╟─bd1fb060-786b-11eb-076a-998aee8fa485
# ╟─76227a6c-7870-11eb-32b3-8392e158059b
# ╟─ca278bd6-89a3-11eb-2388-1d50ae560b7c
# ╠═c1cc1e72-89a6-11eb-3c1d-8ba9aa0b5bb4
# ╠═cbb85948-89a6-11eb-092e-f1194a9774d9
# ╟─f0141896-89a6-11eb-05fe-9d140d242105
# ╟─64af7e68-89b1-11eb-3d43-c1a80cb69dcd
# ╠═6cf7c3e6-89b1-11eb-3374-83471f184496
# ╠═8d9f1b94-89b1-11eb-12c2-991d640e839b
# ╟─e84b96f2-89a3-11eb-3f92-4d4b2849b2e0
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
# ╠═f0b07c72-8999-11eb-316a-293e6cee3c88
# ╟─e39c24c8-8999-11eb-0def-7be2b449dd4f
# ╟─93abda52-786d-11eb-2b2e-0787a202c609
# ╠═b37ec12a-786d-11eb-286c-01ba8c3546c8
# ╟─695e405c-786d-11eb-0a6e-bb776d9626ad
