```@raw html
<style>
    table {
        display: table !important;
        margin: 2rem auto !important;
        border-top: 2pt solid rgba(0,0,0,0.2);
        border-bottom: 2pt solid rgba(0,0,0,0.2);
    }

    pre, div {
        margin-top: 1.4rem !important;
        margin-bottom: 1.4rem !important;
    }

    .code-output {
        padding: 0.7rem 0.5rem !important;
    }

    .admonition-body {
        padding: 0em 1.25em !important;
    }
</style>

<!-- PlutoStaticHTML.Begin -->
<!--
    # This information is used for caching.
    [PlutoStaticHTML.State]
    input_sha = "905613ac95af06caaa111e3f308ddc765738b44fb6f8257b9c40a076d8d173e6"
    julia_version = "1.8.0-rc1"
-->
<pre class='language-julia'><code class='language-julia'>html"""
&lt;div align="center"&gt;
&lt;a class="Header-link " href="https://github.com/TensorBFS/TropicalGEMM.jl" data-hotkey="g d" aria-label="Homepage " data-ga-click="Header, go to dashboard, icon:logo"&gt;
  &lt;svg class="octicon octicon-mark-github v-align-middle" height="32" viewBox="0 0 16 16" version="1.1" width="32" aria-hidden="true"&gt;&lt;path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"&gt;&lt;/path&gt;&lt;/svg&gt;
&lt;/a&gt;
&lt;br&gt;
&lt;a href="https://raw.githubusercontent.com/GiggleLiu/notebooks/master/notebooks/tropical/tropicalgemm.jl" target="_blank"&gt; download this notebook &lt;/a&gt;&lt;/div&gt;
"""</code></pre>
<div align="center">
<a class="Header-link " href="https://github.com/TensorBFS/TropicalGEMM.jl" data-hotkey="g d" aria-label="Homepage " data-ga-click="Header, go to dashboard, icon:logo">
  <svg class="octicon octicon-mark-github v-align-middle" height="32" viewBox="0 0 16 16" version="1.1" width="32" aria-hidden="true"><path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"></path></svg>
</a>
<br>
<a href="https://raw.githubusercontent.com/GiggleLiu/notebooks/master/notebooks/tropical/tropicalgemm.jl" target="_blank"> download this notebook </a></div>



<div class="markdown"><h1>Speed up Tropical matrix multiplication</h1>
</div>


<div class="markdown"><p>By: GiggleLiu and Chris Elrod</p>
</div>


<div class="markdown"><p>This blog is about how to make a GEMM extension for Tropical numbers &#40;<a href="https://github.com/TensorBFS/TropicalGEMM.jl/">TropicalGEMM.jl</a>&#41;, with a close to theoretical optimal performance. It is based on</p>
<ul>
<li><p><a href="https://github.com/JuliaSIMD/LoopVectorization.jl/">LoopVectorization.jl</a>, for vectorizing loops &#40;i.e. utilizing SIMD&#41;,</p>
</li>
<li><p>and <a href="https://github.com/JuliaLinearAlgebra/Octavian.jl">Octavian.jl</a>, a native Julia GEMM library with similar to MKL performance.</p>
</li>
</ul>
</div>


<div class="markdown"><p>Tropical numbers are numbers with tropical algebra. Tropical algebra is defined by replacing the usual sum and product operators for ordinary real numbers with the max and sum operators respectively</p>
<p class="tex">$$\begin&#123;align&#125;
&amp;a ⊕ b &#61; \max&#40;a, b&#41;\\
&amp;a ⊙ b &#61; a &#43; b
\end&#123;align&#125;$$</p>
</div>


<div class="markdown"><p>Its zero and one elements are mapped to regular <code>-Inf</code> and <code>0</code>. For someone who wants to know more about how tropical GEMM can be useful, we highly recommend reading another pluto notebook  <div align=center> <a href='https://giggleliu.github.io/notebooks/tropical/tropicaltensornetwork.html' target=blank>Tropical tensor networks</a> </div> </p>
</div>


<div class="markdown"><p>We choose LoopVectorization and Octavian because it is <strong>fast</strong> and is written in <strong>pure julia</strong>. It has devided the matrix multiplication into small pieces, so that we do not need to handle technical details such as tiling. What people need to do is just implementing several interfaces.</p>
</div>


<div class="markdown"><h2>Let&#39;s jump to the Benchmarks</h2>
</div>


<div class="markdown"><p>The goal is to sqeeze every drop of its computing power of our computing device for testing <code>Intel&#40;R&#41; Core&#40;TM&#41; i5-10400 CPU @ 2.90GHz</code>.  Its theoretical serial computing power for computing a Float64 matrix multiplier is</p>
<pre><code>Serial CPU power &#61; 2.9 GHz &#40;CPU clock speed, we use the maximum Turbo frequency&#41;
			  * 2 &#40;multiplication and add can happen at the same CPU clock&#41;
			  * 2 &#40;number of instructions per cycle&#41;
		      * 4 &#40;avx instruction set has a 256 with register, it can
                   crunch 4 vectorized double precision floating point
				   operations at one CPU cycle&#41;
			&#61; 46.4 GFLOPS</code></pre>
</div>


<div class="markdown"><p>However, the theoretical computing power for tropical matrix multplication is half of that for floating point numbers, because it does not have <code>fma</code> like shortcut to do <code>*</code> and <code>&#43;</code> in a same CPU cycle. So the theoretical maximum computing power for the TropicalF64 GEMM is <code>23.2 GFLOPS</code>.</p>
</div>


<div class="markdown"><p>For matrix size <code>n x n</code>, we show the benchmark results below</p>
<p><img src="https://github.com/TensorBFS/TropicalGEMM.jl/raw/master/benchmarks/benchmark-float64.png" alt="" /></p>
</div>


<div class="markdown"><p>Check the the <a href="https://github.com/TensorBFS/TropicalGEMM.jl/tree/master/benchmarks">benchmarks folder</a> of TropicalGEMM for more benchmarks of different types.</p>
</div>


<div class="markdown"><h2>Implementations</h2>
</div>


<div class="markdown"><p>We are not going to paste the source code and show how it is implemented in detail, because the source code is available in TropicalGEMM <a href="https://github.com/TensorBFS/TropicalGEMM.jl/blob/master/src/gemm.jl">repo</a>. This chapter wants to show some important concepts and the meaning of interfaces that we overwrote. In the following, unless specified, the interfaces overwriten are from <code>VectorizationBase</code>.</p>
</div>


<div class="markdown"><h6 align=center><span style='background-color:yellow'>Warnings before reading</span></h6>
<p>The method introduced to make a BLAS extension is not garanteed to work for other user defined types. The types would have to map 1-1 to native numbers for it to work well, because LoopVectorization assumes that is the case in a way critical to it&#39;s ability to optimize code. So this works for <code>Tropical</code> numbers, but it wouldn&#39;t &#40;for example&#41; <code>Complex</code> or <code>ForwardDiff.Dual</code> numbers, <code>quarternions</code>, or <code>RGB</code> colors. &#40;Chris Elrod: I&#39;ll probably get around to making things like these work eventually using the AbstractInterpreter interface, but the &quot;todo&quot; list before I get there is still quite long.&#41;</p>
</div>


<div class="markdown"><h3>Concepts</h3>
<h5>1. Element types</h5>
<ul>
<li><p><code>Tropical&#123;VectorizationBase.NativeTypes&#125;</code></p>
</li>
<li><p><code>Tropical&#123;&lt;:VectorizationBase.Vec&#125;</code>, a vector of Tropical numbers</p>
</li>
<li><p><code>Tropical&#123;&lt;:VectorizationBase.VecUnroll&#125;</code>, a bundle of <code>Vec</code>s</p>
</li>
</ul>
<p>Here, <code>NativeTypes</code> incldues common floating point numbers, integers, and bit types. Here we use <code>Tropical&#123;&lt;:Vec&#125;</code> to present a vector &#40;that can fit into an SIMD register&#41; of Tropical numbers rather than something like <code>Vec&#123;N, &lt;:Tropical&#125;</code> because <code>Vec</code> is finally processed by SIMD instructions, it can only contain certain <code>NativeTypes</code>.</p>
</div>

<pre class='language-julia'><code class='language-julia'>using LoopVectorization, VectorizationBase, TropicalGEMM</code></pre>


<pre class='language-julia'><code class='language-julia'>vec = Vec(1.0, 2.0, 3.0, 4.0)</code></pre>
<pre id='var-vec' class='code-output documenter-example-output'>VectorizationBase.Vec{4, Float64}<1.0, 2.0, 3.0, 4.0></pre>


<div class="markdown"><p>If you convert it to a tropical number, You will see a small <code>t</code> after it</p>
</div>

<pre class='language-julia'><code class='language-julia'>Tropical(vec)</code></pre>
<pre id='var-hash611660' class='code-output documenter-example-output'>VectorizationBase.Vec{4, Float64}<1.0, 2.0, 3.0, 4.0>ₜ</pre>


<div class="markdown"><p>The same applies for <code>VecUnroll</code></p>
</div>

<pre class='language-julia'><code class='language-julia'>vec_unroll = VecUnroll((vec, vec))</code></pre>
<pre id='var-vec_unroll' class='code-output documenter-example-output'>2 x VectorizationBase.Vec{4, Float64}
VectorizationBase.Vec{4, Float64}<1.0, 2.0, 3.0, 4.0>
VectorizationBase.Vec{4, Float64}<1.0, 2.0, 3.0, 4.0></pre>

<pre class='language-julia'><code class='language-julia'>Tropical(vec_unroll)</code></pre>
<pre id='var-hash127402' class='code-output documenter-example-output'>2 x VectorizationBase.Vec{4, Float64}
VectorizationBase.Vec{4, Float64}<1.0, 2.0, 3.0, 4.0>
VectorizationBase.Vec{4, Float64}<1.0, 2.0, 3.0, 4.0>ₜ</pre>


<div class="markdown"><p><code>VecUnroll</code> is a vectorized <code>Vec</code>. The reason why we need <code>VecUnroll</code> is because it is often faster to unroll a small bundle of vectorized instructions in a loop.</p>
</div>


<div class="markdown"><h5>2. Masks</h5>
<p>A mask is mainly used to avoid loading/storing elements out of bounds &#40;Q: it is correct to say out of bounds?&#41;. When overload an interface, we often implement both the masked and the non-masked versions.</p>
</div>


<div class="markdown"><p>Q: What is the EVLMask?</p>
</div>

<pre class='language-julia'><code class='language-julia'>subtypes(VectorizationBase.AbstractMask)</code></pre>
<pre id='var-hash992618' class='code-output documenter-example-output'>2-element Vector{Any}:
 EVLMask{W, U} where {W, U<:Union{UInt128, UInt16, UInt32, UInt64, UInt8}}
 Mask{W, U} where {W, U<:Union{UInt128, UInt16, UInt32, UInt64, UInt8}}</pre>

<pre class='language-julia'><code class='language-julia'>m = VectorizationBase.Mask{4}(0xe)</code></pre>
<pre id='var-m' class='code-output documenter-example-output'>Mask{4,Bit}<0, 1, 1, 1></pre>


<div class="markdown"><h5>3. Indices</h5>
</div>


<div class="markdown"><p>There are various types of indices</p>
</div>

<pre class='language-julia'><code class='language-julia'>VectorizationBase.Index</code></pre>
<pre id='var-hash592348' class='code-output documenter-example-output'>Union{Int16, Int32, Int64, Int8, UInt16, UInt32, UInt64, UInt8, StaticInt, VectorizationBase.LazyMulAdd{<:Any, <:Any, <:Union{Int16, Int32, Int64, Int8, UInt16, UInt32, UInt64, UInt8, StaticInt}}, Union{VectorizationBase.LazyMulAdd{<:Any, <:Any, <:Union{MM{W}, Unroll{<:Any, <:Any, <:Any, <:Any, W}, Vec{W}}}, MM{W}, Unroll{<:Any, <:Any, <:Any, <:Any, W}, Vec{W}} where W}</pre>


<div class="markdown"><p>One can use the <code>MM</code> type to load a vectorized data into the SIMD register. For example, To continuously load 4 double precision floating point number &#40;8 bytes&#41; from position 0 into a <code>Vec</code>, we can use the following index</p>
</div>

<pre class='language-julia'><code class='language-julia'>vec_index = MM(StaticInt(4), StaticInt(0), StaticInt(8))</code></pre>
<pre id='var-vec_index' class='code-output documenter-example-output'>VectorizationBase.MM{4, 8, Static.StaticInt{0}}<static(0), static(8), static(16), static(24)></pre>


<div class="markdown"><h3>Interfaces to overwrite</h3>
</div>


<div class="markdown"><h4>1. Tell <code>@avx</code> Tropical numbers are is compatible with SIMD</h4>
</div>


<div class="markdown"><ul>
<li><p><code>LoopVectorization.check_args</code> and <code>LoopVectorization.check_type</code></p>
</li>
</ul>
</div>


<div class="markdown"><p>The first thing is telling <code>@avx</code> macro that the <code>Tropical</code> type can utilize SIMD to avoid running into the fallback implementations. <code>@avx</code> is a macro in <code>LoopVectorization</code> that designed to vectorize a loop automatically, it is the corner stone of <code>Octavian</code>.</p>
<p>The <code>@avx</code> macro also checks the array arguments using <code>LoopVectorization.check_args</code> to try and determine if they are compatible with the macro. If <code>check_args</code> returns false, a fall back loop annotated with <code>@inbounds</code> and <code>@fastmath</code> is generated. Note that <code>VectorizationBase</code> provides functions such as <code>vadd</code> and <code>vmul</code> that will ignore <code>@fastmath</code>, preserving IEEE semantics both within <code>@avx</code> and <code>@fastmath</code>. <code>check_args</code> currently returns false for some wrapper types like <code>LinearAlgebra.UpperTriangular</code>, requiring you to use their <code>parent</code>. Triangular loops aren&#39;t yet supported.</p>
<p>LoopVectorization will optimize an <code>@avx</code> loop if <code>check_args</code> on each on the indexed abstract arrays returns true. It returns true for <code>AbstractArray&#123;T&#125;</code>s when <code>check_type&#40;T&#41; &#61;&#61; true</code> and the array or its parent is a <code>StridedArray</code> or <code>AbstractRange</code>.</p>
<p>To provide support for a custom array type, ensure that <code>check_args</code> returns true, either through overloading it or subtyping <code>DenseArray</code>. Additionally, define <code>pointer</code> and <code>stride</code> methods.</p>
</div>

<pre class='language-julia'><code class='language-julia'>LoopVectorization.check_args(TropicalF64, TropicalF64)</code></pre>
<pre id='var-hash109944' class='code-output documenter-example-output'>true</pre>

<pre class='language-julia'><code class='language-julia'>LoopVectorization.check_type(TropicalF64)</code></pre>
<pre id='var-hash436965' class='code-output documenter-example-output'>true</pre>


<div class="markdown"><h4>2. Storage manipulation</h4>
</div>

<pre class='language-julia'><code class='language-julia'>v = Tropical.(randn(10))</code></pre>
<pre id='var-v' class='code-output documenter-example-output'>10-element Vector{TropicalF64}:
  -0.2945481602054801ₜ
   1.0190885428336303ₜ
  -1.4212223893022289ₜ
 -0.11923623603233438ₜ
  -1.9239096962618185ₜ
    0.700790066133385ₜ
  -1.7944285769195838ₜ
 -0.39814213553875527ₜ
   1.4269238867332714ₜ
  -0.8852115298016074ₜ</pre>


<div class="markdown"><ul>
<li><p><code>stridedpointer</code> and <code>gep</code> &#40;create pointers&#41;</p>
</li>
</ul>
</div>


<div class="markdown"><p>e.g. create a strided pointer for an array</p>
</div>

<pre class='language-julia'><code class='language-julia'>ptr = VectorizationBase.stridedpointer(v)</code></pre>
<pre id='var-ptr' class='code-output documenter-example-output'>LayoutPointers.StridedPointer{TropicalF64, 1, 1, 0, (1,), Tuple{StaticInt{8}}, Tuple{StaticInt{1}}}(Ptr{TropicalF64} @0x00007f83b92467a0, ArrayInterface.StrideIndex{1, (1,), 1, Tuple{StaticInt{8}}, Tuple{StaticInt{1}}}((static(8),), (static(1),)))</pre>


<div class="markdown"><p>???</p>
</div>

<pre class='language-julia'><code class='language-julia'>VectorizationBase.gep(ptr.p, 1)</code></pre>
<pre id='var-hash829510' class='code-output documenter-example-output'>Ptr{TropicalF64} @0x00007f83b92467a1</pre>


<div class="markdown"><ul>
<li><p><code>_vload</code> and <code>__vload</code> &#40;loading data&#41;</p>
</li>
</ul>
</div>

<pre class='language-julia'><code class='language-julia'>VectorizationBase.vload(ptr, (3,))</code></pre>
<pre id='var-hash847646' class='code-output documenter-example-output'>-1.4212223893022289ₜ</pre>


<div class="markdown"><p>e.g. load data into a 32*8 bit long register, and the offsets are &#40;0, 8, 16, 24&#41; bits, and mask out the first value.</p>
</div>

<pre class='language-julia'><code class='language-julia'>vi = MM(StaticInt(4), StaticInt(0), StaticInt(8))</code></pre>
<pre id='var-vi' class='code-output documenter-example-output'>VectorizationBase.MM{4, 8, Static.StaticInt{0}}<static(0), static(8), static(16), static(24)></pre>


<div class="markdown"><p>Q: what is 8? why it is used as the normal stride in vload?</p>
</div>

<pre class='language-julia'><code class='language-julia'>VectorizationBase.__vload(ptr.p, vi, m, VectorizationBase.StaticBool(false), StaticInt(32))</code></pre>
<pre id='var-hash173983' class='code-output documenter-example-output'>VectorizationBase.Vec{4, Float64}<0.0, 1.0190885428336303, -1.4212223893022289, -0.11923623603233438>ₜ</pre>


<div class="markdown"><p>If you want to create some zeros</p>
</div>


<div class="markdown"><ul>
<li><p><code>_zero</code> and <code>zero_vecunroll</code> &#40;creating vectorized zero&#41;</p>
</li>
</ul>
</div>


<div class="markdown"><p>e.g. create a vectorized zero of length 4, SMID register size 32 bytes.</p>
</div>

<pre class='language-julia'><code class='language-julia'>VectorizationBase._vzero(StaticInt(4), TropicalF64, StaticInt(32))</code></pre>
<pre id='var-hash341691' class='code-output documenter-example-output'>VectorizationBase.Vec{4, Float64}<-Inf, -Inf, -Inf, -Inf>ₜ</pre>


<div class="markdown"><p>e.g. create 2 vectorized zeros of length 4, SMID register size 32 bytes.</p>
</div>

<pre class='language-julia'><code class='language-julia'> VectorizationBase.zero_vecunroll(StaticInt(2), StaticInt(4), TropicalF64, StaticInt(32))</code></pre>
<pre id='var-hash155212' class='code-output documenter-example-output'>2 x VectorizationBase.Vec{4, Float64}
VectorizationBase.Vec{4, Float64}<-Inf, -Inf, -Inf, -Inf>
VectorizationBase.Vec{4, Float64}<-Inf, -Inf, -Inf, -Inf>ₜ</pre>


<div class="markdown"><ul>
<li><p><code>_vbroadcast</code> &#40;broadcast a scalar to a vector&#41;</p>
</li>
</ul>
</div>


<div class="markdown"><p>e.g. broadcast <code>Tropical&#40;3.0&#41;</code> to SIMD register of size 32 bytes</p>
</div>

<pre class='language-julia'><code class='language-julia'>VectorizationBase._vbroadcast(StaticInt(4), Tropical(3.0), StaticInt(32))</code></pre>
<pre id='var-hash445137' class='code-output documenter-example-output'>VectorizationBase.Vec{4, Float64}<3.0, 3.0, 3.0, 3.0>ₜ</pre>


<div class="markdown"><ul>
<li><p><code>_vstore&#33;</code> and <code>__vstore&#33;</code> &#40;storing data&#41;</p>
</li>
</ul>
</div>


<div class="markdown"><p>e.g. storing a vectorized data into the begining of a vector</p>
</div>

<pre class='language-julia'><code class='language-julia'>let
    v = Tropical.(randn(6))
    ptr = stridedpointer(v)
    vi = MM(StaticInt(4), StaticInt(1), StaticInt(1))
    vstore!(ptr, Tropical(Vec(1.0, 2.0, 3.0, 4.0)), (vi,))
    v
end</code></pre>
<pre id='var-hash756124' class='code-output documenter-example-output'>6-element Vector{TropicalF64}:
                 1.0ₜ
                 2.0ₜ
                 3.0ₜ
                 4.0ₜ
  0.6973826833330837ₜ
 -0.4303812560887033ₜ</pre>


<div class="markdown"><ul>
<li><p><code>Base.promote</code> &#40;promote <code>Tropical&#123;&lt;:Vec&#125;</code> and <code>Tropical&#123;&lt;:VecUnroll&#125;</code>&#41;</p>
</li>
</ul>
</div>


<div class="markdown"><p>e.g. promote <code>Tropical&#40;vec&#41;</code> and <code>Tropical&#40;vec_unroll&#41;</code></p>
</div>

<pre class='language-julia'><code class='language-julia'>promote(Tropical(vec), Tropical(vec_unroll))</code></pre>
<pre id='var-hash175647' class='code-output documenter-example-output'>(2 x VectorizationBase.Vec{4, Float64}
VectorizationBase.Vec{4, Float64}<1.0, 2.0, 3.0, 4.0>
VectorizationBase.Vec{4, Float64}<1.0, 2.0, 3.0, 4.0>ₜ, 2 x VectorizationBase.Vec{4, Float64}
VectorizationBase.Vec{4, Float64}<1.0, 2.0, 3.0, 4.0>
VectorizationBase.Vec{4, Float64}<1.0, 2.0, 3.0, 4.0>ₜ)</pre>


<div class="markdown"><h4>3. Vectorized arithematics</h4>
</div>

<pre class='language-julia'><code class='language-julia'>vec1, vec2, vec3, vec4 = Tropical(Vec(7.0,8.0,3.0,2.0)), Tropical(Vec(1.0,2.0,3.0,4.0)), Tropical(Vec(2.0,2.0,3.0,0.0)), Tropical(Vec(2.0,1.0,1.0,0.0))</code></pre>
<pre id='var-vec4' class='code-output documenter-example-output'>(VectorizationBase.Vec{4, Float64}<7.0, 8.0, 3.0, 2.0>ₜ, VectorizationBase.Vec{4, Float64}<1.0, 2.0, 3.0, 4.0>ₜ, VectorizationBase.Vec{4, Float64}<2.0, 2.0, 3.0, 0.0>ₜ, VectorizationBase.Vec{4, Float64}<2.0, 1.0, 1.0, 0.0>ₜ)</pre>

<pre class='language-julia'><code class='language-julia'>vu = VecUnroll((vec1, vec2, vec3, vec4))</code></pre>
<pre id='var-vu' class='code-output documenter-example-output'>4 x VectorizationBase.Vec{4, Float64}
VectorizationBase.Vec{4, Float64}<7.0, 8.0, 3.0, 2.0>
VectorizationBase.Vec{4, Float64}<1.0, 2.0, 3.0, 4.0>
VectorizationBase.Vec{4, Float64}<2.0, 2.0, 3.0, 0.0>
VectorizationBase.Vec{4, Float64}<2.0, 1.0, 1.0, 0.0>ₜ</pre>


<div class="markdown"><ul>
<li><p><code>Base.FastMath.add_fast</code>, <code>collapse_add</code>, <code>contract_add</code>, <code>reduced_add</code> and <code>vsum</code> &#40;vectorized add&#41;</p>
</li>
</ul>
</div>


<div class="markdown"><p>e.g. <code>vec1 &#43; vec2</code></p>
</div>

<pre class='language-julia'><code class='language-julia'>Base.FastMath.add_fast(vec1, vec2)</code></pre>
<pre id='var-hash122931' class='code-output documenter-example-output'>VectorizationBase.Vec{4, Float64}<7.0, 8.0, 3.0, 4.0>ₜ</pre>


<div class="markdown"><p>We need to handle static integers 0 and 1. They will be used in matrix multiplication as zero and one elements.</p>
</div>

<pre class='language-julia'><code class='language-julia'>Base.FastMath.add_fast(StaticInt(0), vec1)</code></pre>
<pre id='var-hash260928' class='code-output documenter-example-output'>VectorizationBase.Vec{4, Float64}<7.0, 8.0, 3.0, 2.0>ₜ</pre>


<div class="markdown"><p>e.g. <code>&#43;&#40;vu...&#41;</code></p>
</div>

<pre class='language-julia'><code class='language-julia'>VectorizationBase.collapse_add(vu)</code></pre>
<pre id='var-hash256016' class='code-output documenter-example-output'>VectorizationBase.Vec{4, Float64}<7.0, 8.0, 3.0, 4.0>ₜ</pre>


<div class="markdown"><p>e.g. <code>&#40;vec1, vec2, vec3, vec4&#41; -&gt; &#40;vec1&#43;vec2, vec3&#43;vec4&#41;</code></p>
</div>

<pre class='language-julia'><code class='language-julia'>VectorizationBase.contract_add(vu, StaticInt(2))</code></pre>
<pre id='var-hash904083' class='code-output documenter-example-output'>2 x VectorizationBase.Vec{4, Float64}
VectorizationBase.Vec{4, Float64}<7.0, 8.0, 3.0, 2.0>
VectorizationBase.Vec{4, Float64}<2.0, 2.0, 3.0, 4.0>ₜ</pre>


<div class="markdown"><p>e.g. <code>vec1 &#43; vec2</code> &#40;Q: same as add_fast?&#41;</p>
</div>

<pre class='language-julia'><code class='language-julia'>VectorizationBase.reduced_add(vec1, vec2)</code></pre>
<pre id='var-hash123164' class='code-output documenter-example-output'>VectorizationBase.Vec{4, Float64}<7.0, 8.0, 3.0, 4.0>ₜ</pre>


<div class="markdown"><p>e.g. <code>sum&#40;vec1&#41;</code></p>
</div>

<pre class='language-julia'><code class='language-julia'>VectorizationBase.vsum(vec1)</code></pre>
<pre id='var-hash176717' class='code-output documenter-example-output'>8.0ₜ</pre>


<div class="markdown"><ul>
<li><p><code>FastMath.mul_fast</code> &#40;fast multiplication&#41;</p>
</li>
</ul>
</div>


<div class="markdown"><p>e.g. <code>vec1 * vec2</code></p>
</div>

<pre class='language-julia'><code class='language-julia'>Base.FastMath.mul_fast(vec1, vec3)</code></pre>
<pre id='var-hash112191' class='code-output documenter-example-output'>VectorizationBase.Vec{4, Float64}<9.0, 10.0, 6.0, 2.0>ₜ</pre>


<div class="markdown"><p>Handle the one elements properly</p>
</div>

<pre class='language-julia'><code class='language-julia'>Base.FastMath.mul_fast(StaticInt(1), vec1)</code></pre>
<pre id='var-hash322190' class='code-output documenter-example-output'>VectorizationBase.Vec{4, Float64}<7.0, 8.0, 3.0, 2.0>ₜ</pre>


<div class="markdown"><ul>
<li><p><code>Base.fma</code> &#40;fast multiply-add operation&#41;</p>
</li>
</ul>
</div>


<div class="markdown"><p>e.g. Compute <code>vec3 * vec2 &#43; vec1</code></p>
</div>

<pre class='language-julia'><code class='language-julia'>Base.fma(vec3, vec2, vec1)</code></pre>
<pre id='var-hash966795' class='code-output documenter-example-output'>VectorizationBase.Vec{4, Float64}<7.0, 8.0, 6.0, 4.0>ₜ</pre>


<div class="markdown"><h4>4. Other interfaces</h4>
</div>


<div class="markdown"><ul>
<li><p><code>ifelse</code> &#40;vectorized branching&#41;</p>
</li>
</ul>
</div>


<div class="markdown"><p>e.g. <code>masked ? vfmadd_fast&#40;vec1, vec2, vec3&#41; : vec3</code></p>
</div>

<pre class='language-julia'><code class='language-julia'>VectorizationBase.ifelse(VectorizationBase.vfmadd_fast, Mask{4}(0x0e), vec1, vec2, vec3)</code></pre>
<pre id='var-hash169320' class='code-output documenter-example-output'>VectorizationBase.Vec{4, Float64}<2.0, 10.0, 6.0, 6.0>ₜ</pre>


<div class="markdown"><ul>
<li><p><code>vecmaybe</code> &#40;???&#41;</p>
</li>
</ul>
</div>


<div class="markdown"><h2>Comments</h2>
</div>

<pre class='language-julia'><code class='language-julia'>let
    link = html"&lt;div align=center&gt;&lt;a href='https://giggleliu.github.io/notebooks/tropical/tropicaltensornetwork.html' target=blank&gt;Tropical tensor networks&lt;/a&gt;
&lt;/div&gt;"
    md"""
1. Tropical GEMM can be used to find shortest paths, solve combinatoric optimization problems as well as counting solutions. Check $link

2. It is equally important for tropical tensor networks to handle counting tropical algebra
```math
\begin{align}
(x_1, n_1) \odot (x_2,n_2) &= (x_1 + x_2, n_1\cdot n_2)\\
    (x_1, n_1)\oplus (x_2, n_2) &= \begin{cases}
 (x_1\oplus x_2, \, n_1 + n_2 ) & \text{if $x_1 = x_2$} \\
 (x_1\oplus x_2,\, n_1 ) & \text{if $x_1&gt;x_2$} \\
 (x_1\oplus x_2,\, n_2 )& \text{if $x_1 &lt; x_2$}
 \end{cases}.
\end{align}
```

However, composite types are not yet supported in `LoopVectorization`.

3. If you are not sure whether your own type can be accelerated or not, you can catch `Chris Elrod` in the Julia slack channel `#linear-algebra`, he is a smart apple that can answer any question about speeding up a piece of code. If you are interested in discussing Tropical algebra, feel free to ping me (`JinGuo Liu`).
"""
end</code></pre>
<div class="markdown"><ol>
<li><p>Tropical GEMM can be used to find shortest paths, solve combinatoric optimization problems as well as counting solutions. Check <div align=center><a href='https://giggleliu.github.io/notebooks/tropical/tropicaltensornetwork.html' target=blank>Tropical tensor networks</a>
</div></p>
</li>
<li><p>It is equally important for tropical tensor networks to handle counting tropical algebra</p>
</li>
</ol>
<p class="tex">$$\begin&#123;align&#125;
&#40;x_1, n_1&#41; \odot &#40;x_2,n_2&#41; &amp;&#61; &#40;x_1 &#43; x_2, n_1\cdot n_2&#41;\\
    &#40;x_1, n_1&#41;\oplus &#40;x_2, n_2&#41; &amp;&#61; \begin&#123;cases&#125;
 &#40;x_1\oplus x_2, \, n_1 &#43; n_2 &#41; &amp; \text&#123;if &#36;x_1 &#61; x_2&#36;&#125; \\
 &#40;x_1\oplus x_2,\, n_1 &#41; &amp; \text&#123;if &#36;x_1&gt;x_2&#36;&#125; \\
 &#40;x_1\oplus x_2,\, n_2 &#41;&amp; \text&#123;if &#36;x_1 &lt; x_2&#36;&#125;
 \end&#123;cases&#125;.
\end&#123;align&#125;$$</p>
<p>However, composite types are not yet supported in <code>LoopVectorization</code>.</p>
<ol start="3">
<li><p>If you are not sure whether your own type can be accelerated or not, you can catch <code>Chris Elrod</code> in the Julia slack channel <code>#linear-algebra</code>, he is a smart apple that can answer any question about speeding up a piece of code. If you are interested in discussing Tropical algebra, feel free to ping me &#40;<code>JinGuo Liu</code>&#41;.</p>
</li>
</ol>
</div>
<div class='manifest-versions'>
<p>Built with Julia 1.8.0-rc1 and</p>
LoopVectorization 0.12.118<br>
TropicalGEMM 0.1.8<br>
VectorizationBase 0.21.36
</div>

<!-- PlutoStaticHTML.End -->
```
