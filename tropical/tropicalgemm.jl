### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ d3d8702e-8cf8-405a-9b56-45e4153ee265
using LoopVectorization, VectorizationBase, TropicalGEMM

# ╔═╡ 8fe28e17-e716-4a78-a0ba-d70712598a90
html"""
<div align="center">
<a class="Header-link " href="https://github.com/TensorBFS/TropicalGEMM.jl" data-hotkey="g d" aria-label="Homepage " data-ga-click="Header, go to dashboard, icon:logo">
  <svg class="octicon octicon-mark-github v-align-middle" height="32" viewBox="0 0 16 16" version="1.1" width="32" aria-hidden="true"><path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"></path></svg>
</a>
<br>
<a href="https://raw.githubusercontent.com/GiggleLiu/notebooks/master/notebooks/tropical/tropicalgemm.jl" target="_blank"> download this notebook </a></div>
"""

# ╔═╡ 8d24b3aa-7853-11eb-0be4-23088fd5e70a
md"# Speed up Tropical matrix multiplication"

# ╔═╡ d238ceee-f8de-4e14-8947-636c5f879e8c
md"By: GiggleLiu and Chris Elrod"

# ╔═╡ 56082ee0-898f-11eb-13fc-ab4eb456e479
md"This blog is about how to make a GEMM extension for Tropical numbers ([TropicalGEMM.jl](https://github.com/TensorBFS/TropicalGEMM.jl/)), with a close to theoretical optimal performance. It is based on
* [LoopVectorization.jl](https://github.com/JuliaSIMD/LoopVectorization.jl/), for vectorizing loops (i.e. utilizing SIMD),
* and [Octavian.jl](https://github.com/JuliaLinearAlgebra/Octavian.jl), a native Julia GEMM library with similar to MKL performance.
"

# ╔═╡ f0d1bf49-ae16-492d-82a8-b16ea5a54443
md"""
Tropical numbers are numbers with tropical algebra. Tropical algebra is defined by replacing the usual sum and product operators for ordinary real numbers with the max and sum operators respectively
```math
\begin{align}
&a ⊕ b = \max(a, b)\\
&a ⊙ b = a + b
\end{align}
```
"""

# ╔═╡ 31501e08-899d-11eb-0f3c-d95f668c0990
md"""Its zero and one elements are mapped to regular `-Inf` and `0`. For someone who wants to know more about how tropical GEMM can be useful, we highly recommend reading another pluto notebook
$(html"
<div align=center>
<a href='https://giggleliu.github.io/notebooks/tropical/tropicaltensornetwork.html' target=blank>Tropical tensor networks</a>
</div>
")
"""

# ╔═╡ 79f0565c-f187-406b-ba8f-1692fed4773a
md"""
We choose LoopVectorization and Octavian because it is **fast** and is written in **pure julia**. It has devided the matrix multiplication into small pieces, so that we do not need to handle technical details such as tiling. What people need to do is just implementing several interfaces.
"""

# ╔═╡ dba0b4f6-8993-11eb-1822-c993a037dc6b
md"## Let's jump to the Benchmarks"

# ╔═╡ 9bbaefa0-8993-11eb-37a9-854dea2f12dd
md"""
The goal is to sqeeze every drop of its computing power of our computing device for testing `Intel(R) Core(TM) i5-10400 CPU @ 2.90GHz`.  Its theoretical serial computing power for computing a Float64 matrix multiplier is

	Serial CPU power = 2.9 GHz (CPU clock speed, we use the maximum Turbo frequency)
				  * 2 (multiplication and add can happen at the same CPU clock)
				  * 2 (number of instructions per cycle)
			      * 4 (avx instruction set has a 256 with register, it can
                       crunch 4 vectorized double precision floating point
					   operations at one CPU cycle)
				= 46.4 GFLOPS
"""

# ╔═╡ bd1fb060-786b-11eb-076a-998aee8fa485
md"However, the theoretical computing power for tropical matrix multplication is half of that for floating point numbers, because it does not have `fma` like shortcut to do `*` and `+` in a same CPU cycle. So the theoretical maximum computing power for the TropicalF64 GEMM is `23.2 GFLOPS`.
"

# ╔═╡ 28f83f37-200b-4bc1-9cdb-2a461e4262c9
md"""
For matrix size `n x n`, we show the benchmark results below

![](https://github.com/TensorBFS/TropicalGEMM.jl/raw/master/benchmarks/benchmark-float64.png)
"""

# ╔═╡ a3d3557a-ca11-4d14-8cb4-06630bc6badc
md"Check the the [benchmarks folder](https://github.com/TensorBFS/TropicalGEMM.jl/tree/master/benchmarks) of TropicalGEMM for more benchmarks of different types."

# ╔═╡ 76227a6c-7870-11eb-32b3-8392e158059b
md"""
## Implementations
"""

# ╔═╡ f10051e5-eca1-4b36-9c32-479e5e6ca7f7
md"We are not going to paste the source code and show how it is implemented in detail, because the source code is available in TropicalGEMM [repo](https://github.com/TensorBFS/TropicalGEMM.jl/blob/master/src/gemm.jl). This chapter wants to show some important concepts and the meaning of interfaces that we overwrote. In the following, unless specified, the interfaces overwriten are from `VectorizationBase`."

# ╔═╡ def82aee-898e-11eb-3b8d-2325f3709f73
md"""$(HTML("<h6 align=center><span style='background-color:yellow'>Warnings before reading</span></h6>"))

The method introduced to make a BLAS extension is not garanteed to work for other user defined types. The types would have to map 1-1 to native numbers for it to work well, because LoopVectorization assumes that is the case in a way critical to it's ability to optimize code. So this works for `Tropical` numbers, but it wouldn't (for example) `Complex` or `ForwardDiff.Dual` numbers, `quarternions`, or `RGB` colors. (Chris Elrod: I'll probably get around to making things like these work eventually using the AbstractInterpreter interface, but the "todo" list before I get there is still quite long.)
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

# ╔═╡ df211b54-72ca-4d7a-9b17-ac98c97a40ca
md"If you convert it to a tropical number, You will see a small `t` after it"

# ╔═╡ a224ffd8-07cb-4166-9f59-142a63be7634
Tropical(vec)

# ╔═╡ e89c60f5-77ff-413e-bd20-60d6532307e8
md"The same applies for `VecUnroll`"

# ╔═╡ cbb85948-89a6-11eb-092e-f1194a9774d9
vec_unroll = VecUnroll((vec, vec))

# ╔═╡ 4326faed-3167-4cea-a809-19e531fa9247
Tropical(vec_unroll)

# ╔═╡ f0141896-89a6-11eb-05fe-9d140d242105
md"`VecUnroll` is a vectorized `Vec`. The reason why we need `VecUnroll` is because it is often faster to unroll a small bundle of vectorized instructions in a loop."

# ╔═╡ 64af7e68-89b1-11eb-3d43-c1a80cb69dcd
md"""
##### 2. Masks
A mask is mainly used to avoid loading/storing elements out of bounds (Q: it is correct to say out of bounds?).
When overload an interface, we often implement both the masked and the non-masked versions.
"""

# ╔═╡ 6b40af40-31ff-47e5-8422-833cec5a731c
md"Q: What is the EVLMask?"

# ╔═╡ 31072d6e-159b-41a5-aee7-85110b4520a0
subtypes(VectorizationBase.AbstractMask)

# ╔═╡ 6cf7c3e6-89b1-11eb-3374-83471f184496
m = VectorizationBase.Mask{4}(0xe)

# ╔═╡ ae998b13-b18f-4ef6-ac1b-61cbde6ac008
md"##### 3. Indices"

# ╔═╡ 06381f8b-3d6d-4b6f-b5e0-2d9e82206ff7
md"There are various types of indices"

# ╔═╡ a055b531-6aa6-4132-a24b-3a5ccda127c6
VectorizationBase.Index

# ╔═╡ 1afbd730-9742-485e-87d4-bfab8882010e
md"One can use the `MM` type to load a vectorized data into the SIMD register. For example, To continuously load 4 double precision floating point number (8 bytes) from position 0 into a `Vec`, we can use the following index"

# ╔═╡ d0667419-7f6a-4b37-89a4-3b7302e14cec
vec_index = MM(StaticInt(4), StaticInt(0), StaticInt(8))

# ╔═╡ fac58f09-a812-4ea6-9b22-9207b78b6245
md"### Interfaces to overwrite"

# ╔═╡ 82af0af2-786b-11eb-3a49-97519a15a851
md"#### 1. Tell `@avx` Tropical numbers are is compatible with SIMD"

# ╔═╡ d81c61df-1179-462d-97b6-3ed5b4c114b5
md"* `LoopVectorization.check_args` and `LoopVectorization.check_type`"

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
md"#### 2. Storage manipulation
"

# ╔═╡ 2b016b0b-5f22-4347-8c03-fe2c217f73af
v = Tropical.(randn(10))

# ╔═╡ ae06057b-061c-42f3-8c0f-1b62b25d1b45
md"* `stridedpointer` and `gep` (create pointers)"

# ╔═╡ 7bc5bea3-8ec7-458e-89c5-08ad5e9353eb
md"e.g. create a strided pointer for an array"

# ╔═╡ f4f9167b-8b3d-4afe-bb78-475dea40f38c
ptr = VectorizationBase.stridedpointer(v)

# ╔═╡ 2afc083b-2db0-4401-a3d6-e259e5dee09d
md"???"

# ╔═╡ b16fa1f0-466d-4ab1-aa3e-aaa98ec93e86
VectorizationBase.gep(ptr.p, 1)

# ╔═╡ 2f0e159a-6b66-47c4-ac21-350ceeb3a5be
md"""
* `_vload` and `__vload` (loading data)
"""

# ╔═╡ 9a7b58c9-cc12-4292-adb9-5538af4aa3a8
VectorizationBase.vload(ptr, (3,))

# ╔═╡ 9eebd69a-df78-47c2-827b-42691a1d7025
md"e.g. load data into a 32*8 bit long register, and the offsets are (0, 8, 16, 24) bits, and mask out the first value."

# ╔═╡ 7a80faf5-4fd5-4235-a14e-3e6482b57dfb
vi = MM(StaticInt(4), StaticInt(0), StaticInt(8))

# ╔═╡ f349285b-3374-4dde-9b6b-ee39bbb3b951
md"Q: what is 8? why it is used as the normal stride in vload?"

# ╔═╡ 2479b188-b242-487a-b08c-cde2d0d0468a
VectorizationBase.__vload(ptr.p, vi, m, VectorizationBase.StaticBool(false), StaticInt(32))

# ╔═╡ 5a9a6418-2aa7-4b3c-929f-41ee1cb24e09
md"If you want to create some zeros"

# ╔═╡ 49d630cf-3d84-4c0f-a275-40ea04f0fe7c
md"* `_zero` and `zero_vecunroll` (creating vectorized zero)"

# ╔═╡ 6d456655-6541-4d0d-9db4-d5c90a3ec5db
md"e.g. create a vectorized zero of length 4, SMID register size 32 bytes."

# ╔═╡ 7f7e5472-cb7e-47f8-86e2-1f70e848beb0
VectorizationBase._vzero(StaticInt(4), TropicalF64, StaticInt(32))

# ╔═╡ 3138365b-baeb-40d4-90d2-53acc0088754
md"e.g. create 2 vectorized zeros of length 4, SMID register size 32 bytes."

# ╔═╡ 085f476f-7273-4b3a-af7f-2dd39ae4f803
 VectorizationBase.zero_vecunroll(StaticInt(2), StaticInt(4), TropicalF64, StaticInt(32))

# ╔═╡ 97e1b1e2-11b4-4181-a8c2-2322490675d4
md"* `_vbroadcast` (broadcast a scalar to a vector)"

# ╔═╡ ea819ba7-817e-4e60-a69b-ac971bfc4168
md"e.g. broadcast `Tropical(3.0)` to SIMD register of size 32 bytes"

# ╔═╡ 18e0502d-bd5d-44bd-b553-c9515d46a19d
VectorizationBase._vbroadcast(StaticInt(4), Tropical(3.0), StaticInt(32))

# ╔═╡ c345eb4f-3024-4752-a524-3912470b8567
md"""
* `_vstore!` and `__vstore!` (storing data)
"""

# ╔═╡ 82ac4e98-bbc7-4ed5-8c80-b0b3f36523d4
md"e.g. storing a vectorized data into the begining of a vector"

# ╔═╡ 0846196d-45f0-432b-9bf1-3ef8a96c2f46
let
	v = Tropical.(randn(6))
	ptr = stridedpointer(v)
	vi = MM(StaticInt(4), StaticInt(1), StaticInt(1))
	vstore!(ptr, Tropical(Vec(1.0, 2.0, 3.0, 4.0)), (vi,))
	v
end

# ╔═╡ 7a2d76e6-7179-4132-ad4a-8be73a21798d
md"* `Base.promote` (promote `Tropical{<:Vec}` and `Tropical{<:VecUnroll}`)"

# ╔═╡ dd48596e-12b0-48a6-bafc-87e6c3e0b8ff
md"e.g. promote `Tropical(vec)` and `Tropical(vec_unroll)`"

# ╔═╡ 2a8815b6-3969-46e7-8f51-59c4a606e6d5
promote(Tropical(vec), Tropical(vec_unroll))

# ╔═╡ 00581528-78b2-11eb-0904-e9ab3fd376ce
md"#### 3. Vectorized arithematics"

# ╔═╡ 2be3b482-15ac-46f2-af85-11fa662dffb4
vec1, vec2, vec3, vec4 = Tropical(Vec(7.0,8.0,3.0,2.0)), Tropical(Vec(1.0,2.0,3.0,4.0)), Tropical(Vec(2.0,2.0,3.0,0.0)), Tropical(Vec(2.0,1.0,1.0,0.0))

# ╔═╡ 018867bf-d8b9-4292-92a4-6fad9f8b6231
vu = VecUnroll((vec1, vec2, vec3, vec4))

# ╔═╡ fb08d7fe-2561-468c-a529-0519672b90d6
md"* `Base.FastMath.add_fast`, `collapse_add`, `contract_add`, `reduced_add` and `vsum` (vectorized add)"

# ╔═╡ 6601094e-832f-4428-9985-9285bb199884
md"e.g. `vec1 + vec2`"

# ╔═╡ e7b050ba-067b-44bd-bdc0-21ee3524c08f
Base.FastMath.add_fast(vec1, vec2)

# ╔═╡ 31727a71-8a3a-4a3c-954e-18f1c2bc5256
md"We need to handle static integers 0 and 1. They will be used in matrix multiplication as zero and one elements."

# ╔═╡ f1b2cc2d-4fd1-4ac9-91ff-3524266348cf
Base.FastMath.add_fast(StaticInt(0), vec1)

# ╔═╡ 77b2b4e7-e4bb-46c7-9d6e-b047809eee5d
md"e.g. `+(vu...)`"

# ╔═╡ e197c1df-c8d4-4efc-9781-3214d45c7d81
VectorizationBase.collapse_add(vu)

# ╔═╡ 88a7b7ed-f7f7-4e7b-b094-e441a848a136
md"e.g. `(vec1, vec2, vec3, vec4) -> (vec1+vec2, vec3+vec4)`"

# ╔═╡ 972ea87b-f0f1-4eba-8497-897c8b8d7ac1
VectorizationBase.contract_add(vu, StaticInt(2))

# ╔═╡ 39c88348-65b0-449a-a2c5-9ee527fc3e47
md"e.g. `vec1 + vec2` (Q: same as add_fast?)"

# ╔═╡ cdba5628-6856-4ab0-8baa-dba64e85590b
VectorizationBase.reduced_add(vec1, vec2)

# ╔═╡ eb16eb9e-faa4-49cf-91ab-a5b30c7316df
md"e.g. `sum(vec1)`"

# ╔═╡ da1bc10b-4fc7-48b5-83c6-c39c8337b93a
VectorizationBase.vsum(vec1)

# ╔═╡ 98f00ec0-142b-4cd5-a865-d9940482d9fb
md"* `FastMath.mul_fast` (fast multiplication)"

# ╔═╡ 0e2c4923-c00c-47b9-99d6-946c71dc9835
md"e.g. `vec1 * vec2`"

# ╔═╡ 0416901e-5b0f-4fc8-ac34-b4e9ea7e7e9a
Base.FastMath.mul_fast(vec1, vec3)

# ╔═╡ 36814ae9-af88-4a1b-90f1-9d56d2deece1
md"Handle the one elements properly"

# ╔═╡ 8d416951-0f54-4ba9-953c-c33027529583
Base.FastMath.mul_fast(StaticInt(1), vec1)

# ╔═╡ 0880eec7-c3ea-49c1-908f-a615e3bf3f5d
md"* `Base.fma` (fast multiply-add operation)"

# ╔═╡ 89a1b9d6-963d-485b-b623-8ebf08075ed0
md"e.g. Compute `vec3 * vec2 + vec1`"

# ╔═╡ 90037d86-e61e-40d6-a5d9-3c94ee3f2b71
Base.fma(vec3, vec2, vec1)

# ╔═╡ 4e473cf2-78b2-11eb-2fd3-035e416d1650
md"#### 4. Other interfaces"

# ╔═╡ a0ce6c6c-195d-436d-a13c-605ae4727dbf
md"* `ifelse` (vectorized branching)"

# ╔═╡ 79d22d92-af50-4813-b0e2-6e9af834c731
md"e.g. `masked ? vfmadd_fast(vec1, vec2, vec3) : vec3`"

# ╔═╡ ca3e3907-7615-4ec0-88fe-001268221417
VectorizationBase.ifelse(VectorizationBase.vfmadd_fast, Mask{4}(0x0e), vec1, vec2, vec3)

# ╔═╡ e140431a-7637-4f45-9267-9be1527c4779
md"""
* `vecmaybe` (???)
"""

# ╔═╡ 05650690-a5b0-4828-847d-970497caa1a7
md"""
## Comments
"""

# ╔═╡ 695e405c-786d-11eb-0a6e-bb776d9626ad
let
	link = html"<div align=center><a href='https://giggleliu.github.io/notebooks/tropical/tropicaltensornetwork.html' target=blank>Tropical tensor networks</a>
</div>"
	md"""
1. Tropical GEMM can be used to find shortest paths, solve combinatoric optimization problems as well as counting solutions. Check $link

2. It is equally important for tropical tensor networks to handle counting tropical algebra
```math
\begin{align}
(x_1, n_1) \odot (x_2,n_2) &= (x_1 + x_2, n_1\cdot n_2)\\
    (x_1, n_1)\oplus (x_2, n_2) &= \begin{cases}
 (x_1\oplus x_2, \, n_1 + n_2 ) & \text{if $x_1 = x_2$} \\
 (x_1\oplus x_2,\, n_1 ) & \text{if $x_1>x_2$} \\
 (x_1\oplus x_2,\, n_2 )& \text{if $x_1 < x_2$}
 \end{cases}.
\end{align}
```

However, composite types are not yet supported in `LoopVectorization`.

3. If you are not sure whether your own type can be accelerated or not, you can catch `Chris Elrod` in the Julia slack channel `#linear-algebra`, he is a smart apple that can answer any question about speeding up a piece of code. If you are interested in discussing Tropical algebra, feel free to ping me (`JinGuo Liu`).
"""
end

# ╔═╡ Cell order:
# ╟─8fe28e17-e716-4a78-a0ba-d70712598a90
# ╟─8d24b3aa-7853-11eb-0be4-23088fd5e70a
# ╟─d238ceee-f8de-4e14-8947-636c5f879e8c
# ╟─56082ee0-898f-11eb-13fc-ab4eb456e479
# ╟─f0d1bf49-ae16-492d-82a8-b16ea5a54443
# ╟─31501e08-899d-11eb-0f3c-d95f668c0990
# ╟─79f0565c-f187-406b-ba8f-1692fed4773a
# ╟─dba0b4f6-8993-11eb-1822-c993a037dc6b
# ╟─9bbaefa0-8993-11eb-37a9-854dea2f12dd
# ╟─bd1fb060-786b-11eb-076a-998aee8fa485
# ╟─28f83f37-200b-4bc1-9cdb-2a461e4262c9
# ╟─a3d3557a-ca11-4d14-8cb4-06630bc6badc
# ╟─76227a6c-7870-11eb-32b3-8392e158059b
# ╟─f10051e5-eca1-4b36-9c32-479e5e6ca7f7
# ╟─def82aee-898e-11eb-3b8d-2325f3709f73
# ╟─ca278bd6-89a3-11eb-2388-1d50ae560b7c
# ╠═d3d8702e-8cf8-405a-9b56-45e4153ee265
# ╠═c1cc1e72-89a6-11eb-3c1d-8ba9aa0b5bb4
# ╟─df211b54-72ca-4d7a-9b17-ac98c97a40ca
# ╠═a224ffd8-07cb-4166-9f59-142a63be7634
# ╟─e89c60f5-77ff-413e-bd20-60d6532307e8
# ╠═cbb85948-89a6-11eb-092e-f1194a9774d9
# ╠═4326faed-3167-4cea-a809-19e531fa9247
# ╟─f0141896-89a6-11eb-05fe-9d140d242105
# ╟─64af7e68-89b1-11eb-3d43-c1a80cb69dcd
# ╟─6b40af40-31ff-47e5-8422-833cec5a731c
# ╠═31072d6e-159b-41a5-aee7-85110b4520a0
# ╠═6cf7c3e6-89b1-11eb-3374-83471f184496
# ╟─ae998b13-b18f-4ef6-ac1b-61cbde6ac008
# ╟─06381f8b-3d6d-4b6f-b5e0-2d9e82206ff7
# ╠═a055b531-6aa6-4132-a24b-3a5ccda127c6
# ╟─1afbd730-9742-485e-87d4-bfab8882010e
# ╠═d0667419-7f6a-4b37-89a4-3b7302e14cec
# ╟─fac58f09-a812-4ea6-9b22-9207b78b6245
# ╟─82af0af2-786b-11eb-3a49-97519a15a851
# ╟─d81c61df-1179-462d-97b6-3ed5b4c114b5
# ╟─8167bf86-7852-11eb-0201-1996d24d3015
# ╠═af79e3c7-b5e5-4187-b99b-f009d7b945d2
# ╠═a824f58c-6c77-4a63-8e9a-3bf1e9f4c519
# ╟─d18a1a44-d3d8-429c-80a5-c4ed352bfb0d
# ╠═2b016b0b-5f22-4347-8c03-fe2c217f73af
# ╟─ae06057b-061c-42f3-8c0f-1b62b25d1b45
# ╟─7bc5bea3-8ec7-458e-89c5-08ad5e9353eb
# ╠═f4f9167b-8b3d-4afe-bb78-475dea40f38c
# ╟─2afc083b-2db0-4401-a3d6-e259e5dee09d
# ╠═b16fa1f0-466d-4ab1-aa3e-aaa98ec93e86
# ╟─2f0e159a-6b66-47c4-ac21-350ceeb3a5be
# ╠═9a7b58c9-cc12-4292-adb9-5538af4aa3a8
# ╟─9eebd69a-df78-47c2-827b-42691a1d7025
# ╠═7a80faf5-4fd5-4235-a14e-3e6482b57dfb
# ╟─f349285b-3374-4dde-9b6b-ee39bbb3b951
# ╠═2479b188-b242-487a-b08c-cde2d0d0468a
# ╟─5a9a6418-2aa7-4b3c-929f-41ee1cb24e09
# ╟─49d630cf-3d84-4c0f-a275-40ea04f0fe7c
# ╟─6d456655-6541-4d0d-9db4-d5c90a3ec5db
# ╠═7f7e5472-cb7e-47f8-86e2-1f70e848beb0
# ╟─3138365b-baeb-40d4-90d2-53acc0088754
# ╠═085f476f-7273-4b3a-af7f-2dd39ae4f803
# ╟─97e1b1e2-11b4-4181-a8c2-2322490675d4
# ╟─ea819ba7-817e-4e60-a69b-ac971bfc4168
# ╠═18e0502d-bd5d-44bd-b553-c9515d46a19d
# ╟─c345eb4f-3024-4752-a524-3912470b8567
# ╟─82ac4e98-bbc7-4ed5-8c80-b0b3f36523d4
# ╠═0846196d-45f0-432b-9bf1-3ef8a96c2f46
# ╟─7a2d76e6-7179-4132-ad4a-8be73a21798d
# ╟─dd48596e-12b0-48a6-bafc-87e6c3e0b8ff
# ╠═2a8815b6-3969-46e7-8f51-59c4a606e6d5
# ╟─00581528-78b2-11eb-0904-e9ab3fd376ce
# ╠═2be3b482-15ac-46f2-af85-11fa662dffb4
# ╠═018867bf-d8b9-4292-92a4-6fad9f8b6231
# ╟─fb08d7fe-2561-468c-a529-0519672b90d6
# ╟─6601094e-832f-4428-9985-9285bb199884
# ╠═e7b050ba-067b-44bd-bdc0-21ee3524c08f
# ╟─31727a71-8a3a-4a3c-954e-18f1c2bc5256
# ╠═f1b2cc2d-4fd1-4ac9-91ff-3524266348cf
# ╟─77b2b4e7-e4bb-46c7-9d6e-b047809eee5d
# ╠═e197c1df-c8d4-4efc-9781-3214d45c7d81
# ╟─88a7b7ed-f7f7-4e7b-b094-e441a848a136
# ╠═972ea87b-f0f1-4eba-8497-897c8b8d7ac1
# ╟─39c88348-65b0-449a-a2c5-9ee527fc3e47
# ╠═cdba5628-6856-4ab0-8baa-dba64e85590b
# ╟─eb16eb9e-faa4-49cf-91ab-a5b30c7316df
# ╠═da1bc10b-4fc7-48b5-83c6-c39c8337b93a
# ╟─98f00ec0-142b-4cd5-a865-d9940482d9fb
# ╟─0e2c4923-c00c-47b9-99d6-946c71dc9835
# ╠═0416901e-5b0f-4fc8-ac34-b4e9ea7e7e9a
# ╟─36814ae9-af88-4a1b-90f1-9d56d2deece1
# ╠═8d416951-0f54-4ba9-953c-c33027529583
# ╟─0880eec7-c3ea-49c1-908f-a615e3bf3f5d
# ╟─89a1b9d6-963d-485b-b623-8ebf08075ed0
# ╠═90037d86-e61e-40d6-a5d9-3c94ee3f2b71
# ╟─4e473cf2-78b2-11eb-2fd3-035e416d1650
# ╟─a0ce6c6c-195d-436d-a13c-605ae4727dbf
# ╟─79d22d92-af50-4813-b0e2-6e9af834c731
# ╠═ca3e3907-7615-4ec0-88fe-001268221417
# ╟─e140431a-7637-4f45-9267-9be1527c4779
# ╟─05650690-a5b0-4828-847d-970497caa1a7
# ╟─695e405c-786d-11eb-0a6e-bb776d9626ad
