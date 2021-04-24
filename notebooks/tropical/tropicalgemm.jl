### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ c456b902-7959-11eb-03ba-dd14a2cd5758
using PlutoUI, TropicalNumbers

# ╔═╡ d3d8702e-8cf8-405a-9b56-45e4153ee265
using LoopVectorization, VectorizationBase

# ╔═╡ 81d4b4a4-c53a-4c29-9021-07eda458ec54
using TropicalGEMM

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

# ╔═╡ 8d24b3aa-7853-11eb-0be4-23088fd5e70a
md"# Tropical GEMM"

# ╔═╡ 56082ee0-898f-11eb-13fc-ab4eb456e479
md"This blog is about how to make a GEMM extension for Tropical numbers, with a close to theoretical optimal performance. It is based on
* `LoopVectorization`, for vectorizing loops (i.e. utilizing SIMD),
* and `Octavian`, a native Julia GEMM library with similar to MKL performance.
"

# ╔═╡ 31501e08-899d-11eb-0f3c-d95f668c0990
md"To know more about how tropical GEMM can be useful, we highly recommend readers to read this pluto notebook about [Tropical tensor networks](https://giggleliu.github.io/notebooks/tropical/tropicaltensornetwork.html)."

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
				  * 2 (number of instructions per cycle)
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

# ╔═╡ 28f83f37-200b-4bc1-9cdb-2a461e4262c9
md"""## Benchmarks
Matrix size `n x n`, CPU Intel(R) Core(TM) i5-10400 CPU @ 2.90GHz. Check the the benchmarks folder of TropicalGEMM for more benchmarks of different types.

![](https://github.com/TensorBFS/TropicalGEMM.jl/raw/master/benchmarks/benchmark-float64.png)
"""

# ╔═╡ 76227a6c-7870-11eb-32b3-8392e158059b
md"""
## Implementations
### Why LoopVectorization and Octavian?

It is fast and extensible. It has devided the problems into small pieces, and handled technical details of tiling et al. What people need to do is determining how to load and store data, how to do elementary arithmetic operations.
"""

# ╔═╡ ca278bd6-89a3-11eb-2388-1d50ae560b7c
md"""### Concepts
##### 1. Element types
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

# ╔═╡ ae998b13-b18f-4ef6-ac1b-61cbde6ac008
md"##### 3. Indices"

# ╔═╡ 1afbd730-9742-485e-87d4-bfab8882010e
md"An index that continuously loading 4 double precision floating point number (8 bytes) from position 0."

# ╔═╡ d0667419-7f6a-4b37-89a4-3b7302e14cec
vec_index = MM(StaticInt(4), StaticInt(0), StaticInt(8))

# ╔═╡ e84b96f2-89a3-11eb-3f92-4d4b2849b2e0
md"""
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

# ╔═╡ af79e3c7-b5e5-4187-b99b-f009d7b945d2
LoopVectorization.check_args(TropicalF64, TropicalF64)

# ╔═╡ a824f58c-6c77-4a63-8e9a-3bf1e9f4c519
LoopVectorization.check_type(TropicalF64)

# ╔═╡ d18a1a44-d3d8-429c-80a5-c4ed352bfb0d
md"## Store and load data
"

# ╔═╡ ae06057b-061c-42f3-8c0f-1b62b25d1b45
md"`stridedpointer` and `gep`"

# ╔═╡ 2b016b0b-5f22-4347-8c03-fe2c217f73af
v = Tropical.(randn(10))

# ╔═╡ 7bc5bea3-8ec7-458e-89c5-08ad5e9353eb
md"create a strided pointer and load the 3rd value"

# ╔═╡ f4f9167b-8b3d-4afe-bb78-475dea40f38c
ptr = VectorizationBase.stridedpointer(v)

# ╔═╡ 2afc083b-2db0-4401-a3d6-e259e5dee09d
md"???"

# ╔═╡ b16fa1f0-466d-4ab1-aa3e-aaa98ec93e86
VectorizationBase.gep(ptr.p, 1)

# ╔═╡ 2f0e159a-6b66-47c4-ac21-350ceeb3a5be
md"""
`_vload` and `__vload`
"""

# ╔═╡ 9a7b58c9-cc12-4292-adb9-5538af4aa3a8
VectorizationBase.vload(ptr, (3,))

# ╔═╡ 9eebd69a-df78-47c2-827b-42691a1d7025
md"load data into a 32*8 bit long register, and the offsets are (0, 8, 16, 24) bits, and mask out the first value."

# ╔═╡ 7a80faf5-4fd5-4235-a14e-3e6482b57dfb
vi = MM(StaticInt(4), StaticInt(0), StaticInt(8))

# ╔═╡ 2479b188-b242-487a-b08c-cde2d0d0468a
VectorizationBase.__vload(ptr.p, vi, m, VectorizationBase.StaticBool(false), StaticInt(32))

# ╔═╡ 5a9a6418-2aa7-4b3c-929f-41ee1cb24e09
md"If you want to create some zeros"

# ╔═╡ 49d630cf-3d84-4c0f-a275-40ea04f0fe7c
md"`_zero` and `zero_vecunroll`"

# ╔═╡ 7f7e5472-cb7e-47f8-86e2-1f70e848beb0
VectorizationBase._vzero(StaticInt(4), TropicalF64, StaticInt(32))

# ╔═╡ 085f476f-7273-4b3a-af7f-2dd39ae4f803
 VectorizationBase.zero_vecunroll(StaticInt(2), StaticInt(4), TropicalF64, StaticInt(32))

# ╔═╡ 97e1b1e2-11b4-4181-a8c2-2322490675d4
md"`_vbroadcast`"

# ╔═╡ 18e0502d-bd5d-44bd-b553-c9515d46a19d
VectorizationBase._vbroadcast(StaticInt(4), Tropical(3.0), StaticInt(32))

# ╔═╡ c345eb4f-3024-4752-a524-3912470b8567
md"""
`_vstore!` and `__vstore!`
"""

# ╔═╡ 2a8815b6-3969-46e7-8f51-59c4a606e6d5
promote(vec, vec_unroll)

# ╔═╡ 00581528-78b2-11eb-0904-e9ab3fd376ce
md"#### vectorized operations"

# ╔═╡ 2be3b482-15ac-46f2-af85-11fa662dffb4
vec1, vec2, vec3, vec4 = Tropical(Vec(7.0,8.0,3.0,2.0)), Tropical(Vec(1.0,2.0,3.0,4.0)), Tropical(Vec(2.0,2.0,3.0,0.0)), Tropical(Vec(2.0,1.0,1.0,0.0))

# ╔═╡ 018867bf-d8b9-4292-92a4-6fad9f8b6231
vu = VecUnroll((vec1, vec2, vec3, vec4))

# ╔═╡ e7b050ba-067b-44bd-bdc0-21ee3524c08f
Base.FastMath.add_fast(vec1, vec2)

# ╔═╡ e197c1df-c8d4-4efc-9781-3214d45c7d81
VectorizationBase.collapse_add(vu)

# ╔═╡ 972ea87b-f0f1-4eba-8497-897c8b8d7ac1
VectorizationBase.contract_add(vu, StaticInt(2))

# ╔═╡ cdba5628-6856-4ab0-8baa-dba64e85590b
VectorizationBase.reduced_add(vec1, vec2)

# ╔═╡ 695c8b2f-a795-4952-97bb-50d1163ae880
Mask{4}(0x0e)

# ╔═╡ 89a1b9d6-963d-485b-b623-8ebf08075ed0
md"Compute `vec3 * vec2 + vec1`"

# ╔═╡ 90037d86-e61e-40d6-a5d9-3c94ee3f2b71
Base.fma(vec3, vec2, vec1)

# ╔═╡ 36814ae9-af88-4a1b-90f1-9d56d2deece1
md"Handle zero elements and one elements properly"

# ╔═╡ f1b2cc2d-4fd1-4ac9-91ff-3524266348cf
Base.FastMath.add_fast(StaticInt(0), vec1)

# ╔═╡ 8d416951-0f54-4ba9-953c-c33027529583
Base.FastMath.mul_fast(StaticInt(1), vec1)

# ╔═╡ ca3e3907-7615-4ec0-88fe-001268221417
VectorizationBase.ifelse(VectorizationBase.vfmadd_fast, Mask{4}(0x0e), vec1, vec2, vec3)

# ╔═╡ da1bc10b-4fc7-48b5-83c6-c39c8337b93a
VectorizationBase.vsum(vec1)

# ╔═╡ 4e473cf2-78b2-11eb-2fd3-035e416d1650
md"## Other interfaces"

# ╔═╡ e140431a-7637-4f45-9267-9be1527c4779
md"""
* `vecmaybe`
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
# ╟─8d24b3aa-7853-11eb-0be4-23088fd5e70a
# ╟─56082ee0-898f-11eb-13fc-ab4eb456e479
# ╟─31501e08-899d-11eb-0f3c-d95f668c0990
# ╟─dba0b4f6-8993-11eb-1822-c993a037dc6b
# ╟─9bbaefa0-8993-11eb-37a9-854dea2f12dd
# ╟─bd1fb060-786b-11eb-076a-998aee8fa485
# ╟─28f83f37-200b-4bc1-9cdb-2a461e4262c9
# ╟─76227a6c-7870-11eb-32b3-8392e158059b
# ╟─def82aee-898e-11eb-3b8d-2325f3709f73
# ╟─ca278bd6-89a3-11eb-2388-1d50ae560b7c
# ╠═c1cc1e72-89a6-11eb-3c1d-8ba9aa0b5bb4
# ╠═cbb85948-89a6-11eb-092e-f1194a9774d9
# ╟─f0141896-89a6-11eb-05fe-9d140d242105
# ╟─64af7e68-89b1-11eb-3d43-c1a80cb69dcd
# ╠═6cf7c3e6-89b1-11eb-3374-83471f184496
# ╟─ae998b13-b18f-4ef6-ac1b-61cbde6ac008
# ╟─1afbd730-9742-485e-87d4-bfab8882010e
# ╠═d0667419-7f6a-4b37-89a4-3b7302e14cec
# ╟─e84b96f2-89a3-11eb-3f92-4d4b2849b2e0
# ╟─82af0af2-786b-11eb-3a49-97519a15a851
# ╟─8167bf86-7852-11eb-0201-1996d24d3015
# ╠═d3d8702e-8cf8-405a-9b56-45e4153ee265
# ╠═81d4b4a4-c53a-4c29-9021-07eda458ec54
# ╠═af79e3c7-b5e5-4187-b99b-f009d7b945d2
# ╠═a824f58c-6c77-4a63-8e9a-3bf1e9f4c519
# ╟─d18a1a44-d3d8-429c-80a5-c4ed352bfb0d
# ╟─ae06057b-061c-42f3-8c0f-1b62b25d1b45
# ╠═2b016b0b-5f22-4347-8c03-fe2c217f73af
# ╟─7bc5bea3-8ec7-458e-89c5-08ad5e9353eb
# ╠═f4f9167b-8b3d-4afe-bb78-475dea40f38c
# ╟─2afc083b-2db0-4401-a3d6-e259e5dee09d
# ╠═b16fa1f0-466d-4ab1-aa3e-aaa98ec93e86
# ╟─2f0e159a-6b66-47c4-ac21-350ceeb3a5be
# ╠═9a7b58c9-cc12-4292-adb9-5538af4aa3a8
# ╟─9eebd69a-df78-47c2-827b-42691a1d7025
# ╠═7a80faf5-4fd5-4235-a14e-3e6482b57dfb
# ╠═2479b188-b242-487a-b08c-cde2d0d0468a
# ╟─5a9a6418-2aa7-4b3c-929f-41ee1cb24e09
# ╟─49d630cf-3d84-4c0f-a275-40ea04f0fe7c
# ╠═7f7e5472-cb7e-47f8-86e2-1f70e848beb0
# ╠═085f476f-7273-4b3a-af7f-2dd39ae4f803
# ╟─97e1b1e2-11b4-4181-a8c2-2322490675d4
# ╠═18e0502d-bd5d-44bd-b553-c9515d46a19d
# ╟─c345eb4f-3024-4752-a524-3912470b8567
# ╠═2a8815b6-3969-46e7-8f51-59c4a606e6d5
# ╟─00581528-78b2-11eb-0904-e9ab3fd376ce
# ╠═2be3b482-15ac-46f2-af85-11fa662dffb4
# ╠═018867bf-d8b9-4292-92a4-6fad9f8b6231
# ╠═e7b050ba-067b-44bd-bdc0-21ee3524c08f
# ╠═e197c1df-c8d4-4efc-9781-3214d45c7d81
# ╠═972ea87b-f0f1-4eba-8497-897c8b8d7ac1
# ╠═cdba5628-6856-4ab0-8baa-dba64e85590b
# ╠═695c8b2f-a795-4952-97bb-50d1163ae880
# ╟─89a1b9d6-963d-485b-b623-8ebf08075ed0
# ╠═90037d86-e61e-40d6-a5d9-3c94ee3f2b71
# ╟─36814ae9-af88-4a1b-90f1-9d56d2deece1
# ╠═f1b2cc2d-4fd1-4ac9-91ff-3524266348cf
# ╠═8d416951-0f54-4ba9-953c-c33027529583
# ╠═ca3e3907-7615-4ec0-88fe-001268221417
# ╠═da1bc10b-4fc7-48b5-83c6-c39c8337b93a
# ╟─4e473cf2-78b2-11eb-2fd3-035e416d1650
# ╟─e140431a-7637-4f45-9267-9be1527c4779
# ╟─93abda52-786d-11eb-2b2e-0787a202c609
# ╟─b37ec12a-786d-11eb-286c-01ba8c3546c8
# ╟─695e405c-786d-11eb-0a6e-bb776d9626ad
