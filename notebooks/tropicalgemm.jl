### A Pluto.jl notebook ###
# v0.19.8

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
<a href="https://raw.githubusercontent.com/GiggleLiu/notebooks/master/notebooks/tropicalgemm.jl" target="_blank"> download this notebook </a></div>
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

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LoopVectorization = "bdcacae8-1622-11e9-2a5c-532679323890"
TropicalGEMM = "a4ad3063-64a7-4bad-8738-34ed09bc0236"
VectorizationBase = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"

[compat]
LoopVectorization = "~0.12.118"
TropicalGEMM = "~0.1.8"
VectorizationBase = "~0.21.36"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.0-rc1"
manifest_format = "2.0"
project_hash = "bd0478889d00f12ced3b8fb9934327a69276d741"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["ArrayInterfaceCore", "Compat", "IfElse", "LinearAlgebra", "Static"]
git-tree-sha1 = "dce69568d03b16af8ed09896f27768d8d983d819"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "6.0.15"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "d0f59ebfe8d3ea2799fb3fb88742d69978e5843e"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.10"

[[deps.ArrayInterfaceOffsetArrays]]
deps = ["ArrayInterface", "OffsetArrays", "Static"]
git-tree-sha1 = "7dce0e2846e7496622f5d2742502d7e029693458"
uuid = "015c0d05-e682-4f19-8f0a-679ce4c54826"
version = "0.1.5"

[[deps.ArrayInterfaceStaticArrays]]
deps = ["Adapt", "ArrayInterface", "LinearAlgebra", "Static", "StaticArrays"]
git-tree-sha1 = "d7dc30474e73173a990eca86af76cae8790fa9f2"
uuid = "b0d46f97-bff5-4637-a19a-dd75974142cd"
version = "0.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitTwiddlingConvenienceFunctions]]
deps = ["Static"]
git-tree-sha1 = "28bbdbf0354959db89358d1d79d421ff31ef0b5e"
uuid = "62783981-4cbd-42fc-bca8-16325de8dc4b"
version = "0.1.3"

[[deps.CPUSummary]]
deps = ["CpuId", "IfElse", "Static"]
git-tree-sha1 = "0eaf4aedad5ccc3e39481db55d72973f856dc564"
uuid = "2a0fbf3d-bb9c-48f3-b0a9-814d99fd7ab9"
version = "0.1.22"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9489214b993cd42d17f44c36e359bf6a7c919abf"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.0"

[[deps.CloseOpenIntervals]]
deps = ["ArrayInterface", "Static"]
git-tree-sha1 = "eb61d6b97041496058245821e3bb7eba2b2cf4db"
uuid = "fb6a15b2-703c-40df-9091-08a04967cfa9"
version = "0.1.8"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "924cdca592bc16f14d2f7006754a621735280b74"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.1.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.CpuId]]
deps = ["Markdown"]
git-tree-sha1 = "fcbb72b032692610bfbdb15018ac16a36cf2e406"
uuid = "adafc99b-e345-5852-983c-f28acb93d879"
version = "0.3.1"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[deps.DiffRules]]
deps = ["NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "7220bc21c33e990c14f4a9a319b1d242ebc5b269"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.3.1"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "c5544d8abb854e306b7b2f799ab31cdba527ccae"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "ef3fec65f9db26fa2cf8f4133c697c5b7ce63c1d"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.22"

[[deps.HostCPUFeatures]]
deps = ["BitTwiddlingConvenienceFunctions", "IfElse", "Libdl", "Static"]
git-tree-sha1 = "18be5268cf415b5e27f34980ed25a7d34261aa83"
uuid = "3e5b6fbb-0976-4d2c-9146-d79de83f2fb0"
version = "0.1.7"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.LayoutPointers]]
deps = ["ArrayInterface", "ArrayInterfaceOffsetArrays", "ArrayInterfaceStaticArrays", "LinearAlgebra", "ManualMemory", "SIMDTypes", "Static"]
git-tree-sha1 = "a575de5a424a395217930fea6d0934ea853d0158"
uuid = "10f19ff3-798f-405d-979b-55457f8fc047"
version = "0.1.9"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.81.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoopVectorization]]
deps = ["ArrayInterface", "ArrayInterfaceCore", "ArrayInterfaceOffsetArrays", "ArrayInterfaceStaticArrays", "CPUSummary", "ChainRulesCore", "CloseOpenIntervals", "DocStringExtensions", "ForwardDiff", "HostCPUFeatures", "IfElse", "LayoutPointers", "LinearAlgebra", "OffsetArrays", "PolyesterWeave", "SIMDDualNumbers", "SIMDTypes", "SLEEFPirates", "SpecialFunctions", "Static", "ThreadingUtilities", "UnPack", "VectorizationBase"]
git-tree-sha1 = "5ea9a0aaf5ded7f0b6e43c96ca1793e60c96af93"
uuid = "bdcacae8-1622-11e9-2a5c-532679323890"
version = "0.12.118"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.ManualMemory]]
git-tree-sha1 = "bcaef4fc7a0cfe2cba636d84cda54b5e4e4ca3cd"
uuid = "d125e4d3-2237-4719-b19c-fa641b8a4667"
version = "0.1.8"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Octavian]]
deps = ["ArrayInterface", "CPUSummary", "IfElse", "LoopVectorization", "ManualMemory", "PolyesterWeave", "Requires", "Static", "ThreadingUtilities", "VectorizationBase"]
git-tree-sha1 = "e4705f509d1d623518ac02fdc734e93df980a9df"
uuid = "6fd5a793-0b7e-452c-907f-f8bfe9c57db4"
version = "0.3.14"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "b4975062de00106132d0b01b5962c09f7db7d880"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.5"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PolyesterWeave]]
deps = ["BitTwiddlingConvenienceFunctions", "CPUSummary", "IfElse", "Static", "ThreadingUtilities"]
git-tree-sha1 = "7e597df97e46ffb1c8adbaddfa56908a7a20194b"
uuid = "1d0040c9-8b98-4ee7-8388-3f51789ca0ad"
version = "0.1.5"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMDDualNumbers]]
deps = ["ForwardDiff", "IfElse", "SLEEFPirates", "VectorizationBase"]
git-tree-sha1 = "dd4195d308df24f33fb10dde7c22103ba88887fa"
uuid = "3cdde19b-5bb0-4aaf-8931-af3e248e098b"
version = "0.1.1"

[[deps.SIMDTypes]]
git-tree-sha1 = "330289636fb8107c5f32088d2741e9fd7a061a5c"
uuid = "94e857df-77ce-4151-89e5-788b33177be4"
version = "0.1.0"

[[deps.SLEEFPirates]]
deps = ["IfElse", "Static", "VectorizationBase"]
git-tree-sha1 = "ac399b5b163b9140f9c310dfe9e9aaa225617ff6"
uuid = "476501e8-09a2-5ece-8869-fb82de89a1fa"
version = "0.6.32"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["OpenSpecFun_jll"]
git-tree-sha1 = "7286f31f27e3335cba31c618ac344a35eceac060"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.1.0"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "5d2c08cef80c7a3a8ba9ca023031a85c263012c5"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.6.6"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "2bbd9f2e40afd197a1379aef05e0d85dba649951"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.7"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadingUtilities]]
deps = ["ManualMemory"]
git-tree-sha1 = "f8629df51cab659d70d2e5618a430b4d3f37f2c3"
uuid = "8290d209-cae3-49c0-8002-c8c24d57dab5"
version = "0.5.0"

[[deps.TropicalGEMM]]
deps = ["LinearAlgebra", "LoopVectorization", "Octavian", "TropicalNumbers", "VectorizationBase"]
git-tree-sha1 = "63838f5be6c1f591e124a6df4df55c5ae520a25f"
uuid = "a4ad3063-64a7-4bad-8738-34ed09bc0236"
version = "0.1.8"

[[deps.TropicalNumbers]]
git-tree-sha1 = "f3659ba817a2fdc8665e9bf032d66f1a107a56af"
uuid = "b3a74e9c-7526-4576-a4eb-79c0d4c32334"
version = "0.5.3"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.VectorizationBase]]
deps = ["ArrayInterface", "CPUSummary", "HostCPUFeatures", "IfElse", "LayoutPointers", "Libdl", "LinearAlgebra", "SIMDTypes", "Static"]
git-tree-sha1 = "7d3de169cd221392082a5abc7f363726e1a30628"
uuid = "3d5dd08c-fd9d-11e8-17fa-ed2836048c2f"
version = "0.21.36"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.41.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

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
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
