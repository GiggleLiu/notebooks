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
begin
	using Revise, PlutoUI, CoordinateTransformations, StaticArrays, Rotations, Viznet, Compose
	
	# cubic geometry
	struct Cubic{T}
		grid::NTuple{3,T}
	end
	function Base.getindex(c, i, j, k)
		(i,j,k) .* c.grid
	end
	
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
end;

# ╔═╡ 5bb40ad6-7b33-11eb-0b31-63d5e47fa0e7
using TropicalNumbers,  # tropical number type
		LightGraphs,	# graph operations
		Random,
    	SimpleTensorNetworks  # tensor network contraction

# ╔═╡ 7fa86de2-7ac5-11eb-3d55-c9b6cb57c5d6
html"""
<style>
ul li { margin-bottom: 0.5em; }
body {font-size: 18pt;}
</style>
<script>
document.body.onkeyup = function(e) {
if (e.ctrlKey && e.altKey && e.which == 80) {
    present();
} else if (e.ctrlKey && e.which == 37) {
	var prev_button = document.querySelector(".changeslide.prev");
	prev_button.dispatchEvent(new Event('click'));
} else if (e.ctrlKey && e.which == 39) {
	var prev_button = document.querySelector(".changeslide.next");
	prev_button.dispatchEvent(new Event('click'));
  }
};
document.body.onclick = function(e) {
	if (e.target.tagName == 'BODY'){
		e.preventDefault();
		var prev_button = document.querySelector(".changeslide.next");
		prev_button.dispatchEvent(new Event('click'));
} else if (e.target.tagName == 'PLUTO-SHOULDER'){
	e.preventDefault();
	var prev_button = document.querySelector(".changeslide.prev");
	prev_button.dispatchEvent(new Event('click'));
	}
};
</script>
"""


# ╔═╡ 121b4926-7aba-11eb-30e1-7b8edd4f0166
md"""# Tropical tensor networks for solving combinatoric optimization problems


$(HTML("<br><p><big><strong>Tropical tensor network for ground states of spin glasses</strong></big></p>
<p>Phys. Rev. Lett. (26 January 2021)</p>
<p>Jin-Guo Liu, Lei Wang, and Pan Zhang</p>"))


[arxiv 2008.06888](https://arxiv.org/abs/2008.06888)
"""

# ╔═╡ 265649e4-7abb-11eb-1a02-a3c101cef89d
html"""<h1> Goals</h1>
<ul style="list-style-type: square;">
<li>What is a tropical tensor network?</li>
<li>With tropical tensor networks, solving
<ul>
	<li>spin glass,</li>
	<li>2-satistiability,</li>
	<li>Potts model,</li>
	<li>maximum independent set</li>
</ul>and their counting problems .</li>
</ul>
"""

# ╔═╡ 7bdf517e-79ff-11eb-38a3-49c02d94d943
md"## The Song Shan Lake Spring School (SSSS) Challenge"

# ╔═╡ 88e14ef2-7af1-11eb-23d6-b34b1eff8f87
md"""
In the Buckyball structure as shown in the figure, we attach an ising spin ``s_i=\pm 1`` on each vertex. The neighboring spins interact with an anti-ferromagnetic coupling of unit strength.

1. Get ``\ln Z/N``, where ``N`` is the number of vertices, and
```math
Z = \sum_{\{s\}}e^{-\sum_{i,j \in E}s_i s_j}
```

2. Count the ground state degeneracy.
"""

# ╔═╡ a387b018-79ff-11eb-0383-1fcf82853afc
md"""
We get more than 5 solutions for each quiz, they are available in Github repository [QuantumBFS/SSSS](https://github.com/QuantumBFS/SSSS/blob/master/Challenge.md)
"""

# ╔═╡ 3205a536-7a17-11eb-3473-b71305c96ca4
md"# Today, we have one more solution: Tropical tensor networks"

# ╔═╡ be76e52a-7852-11eb-179b-afbc6efcab55
md"## Tropical algebra"

# ╔═╡ d0b54b76-7852-11eb-2398-0911380fa090
md"""

```math
\begin{align}
&a ⊕ b = \max(a, b)\\
&a ⊙ b = a + b
\end{align}
```
"""

# ╔═╡ 9dbbcb08-7d40-11eb-3d58-d1ccf9a9d803
md"# Can we still define matrix multiplication on Tropical numbers?"

# ╔═╡ 90ee0580-7d38-11eb-36d3-299cefec4ddf
function comment(str)
	md"Sheng Tao: $str"
end;

# ╔═╡ 908952d0-7d0d-11eb-2cf7-6fa6d1fbee88
let
	fail = HTML("""<div style='background-color:transparent'><del>$(md"For each a in ``R`` there exists in R such that ``a ⊕ (−a) = 0``   (that is, ``−a`` is the additive inverse of ``a``).  (not true for tropical numbers)" |> html)</del></div>""")
	md"""
## Yes, Because Tropical algebra is a A Semi-Ring over ``\mathbb R \cup \{-\infty\}``
Ring is a set ``R`` equipped with two binary operations ⊕ (addition) and ⊙ (multiplication)

* is an abebian group under ⊕
    * associative
    * commutative
    * zero element (``-infty``), since  $-\infty \oplus x = x  $ and $-\infty \odot x = -\infty$.
    * $(fail)
* is a monoid under ⊙
    * associative
    * There is an element ``1`` in ``R`` such that ``a ⊙ 1 = a`` and ``1 ⊙ a = a``. For tropical numbers, $0$ acts as the multiplicative identity since $0 \odot x = x$.
* ⊙ is distributive with respect to ⊕,
"""
end

# ╔═╡ 98ae0960-797d-11eb-3646-c5b7e05d3f7c
md"""
### Example: Tropical δ tensor
The tropical $\delta$ tensor of rank $n$ and dimension $q$ is defined as
```math
δ_{s_i s_j\ldots s_n}^{n,q} = \begin{cases}
 0, & s_i = s_j =\ldots s_n\\
 -\infty, &otherwise
\end{cases}
```
where $s_i,s_j,\ldots s_n \in \{1,2,\ldots q\}$.
"""

# ╔═╡ 86921d00-7a17-11eb-2695-add5f9eeda5b
md"## Tropical matrix multiplication
```math
\begin{align}
\cancel{Y[i,j] := \sum_k A[i,k] \times B[k,j]}\\

Y[i,j] := \max_k (A[i,k] + B[k,j])
\end{align}
```
"

# ╔═╡ 915e8096-7a17-11eb-177d-a39ffed7ca91
let
	Compose.set_default_graphic_size(15cm, 10/3*cm)
	sq = nodestyle(:circle; r=0.1)
	eb = bondstyle(:line)
	tb = textstyle(:default, fontsize(25px))
	tb2 = textstyle(:default, fontsize(40px), fill("white"))
	y0 = 0.15
	x = (0.3, y0)
	y = (0.7, y0)
	img = canvas() do
		sq >> x
		sq >> y
		eb >> (x, y)
		eb >> (x, (0.0, y0))
		eb >> (x, (1.0, y0))
		tb >> ((0.1, y0+0.05), "i")
		tb >> ((0.9, y0+0.05), "j")
		tb >> ((0.5, y0+0.05), "k")
		tb2 >> (x, "A")
		tb2 >> (y, "B")
	end
	Compose.compose(context(0.38, 0.0, 1/1.5^2, 2.0), img)
end


# ╔═╡ 1cbb86e4-7d39-11eb-05d9-c3ab43a03298
comment("to many words")

# ╔═╡ 958c489e-7a17-11eb-2be4-cfdb44da7d2f
md"
* we use **circles** to represent tropical tensors, lines (or legs) to represent their dimensions,
* connecting a pair of legs means **max-sum** over a specific degree of freedom corresponding to the label,
"

# ╔═╡ ddca618c-7d22-11eb-0c19-e1a149e19f3c
md"#### References
* The complexity of tropical matrix factorization [arxiv:1205.7079](https://arxiv.org/abs/1205.7079)"

# ╔═╡ 3d7ca08c-7b01-11eb-1d78-af35dc7e577c
md"## Example: Tropical matrix multiplication to find the shortest path"

# ╔═╡ 4fbec178-7d39-11eb-1f27-13ac45f53e45
comment("remove")

# ╔═╡ 742fdf74-7d0e-11eb-026f-81ee2c7c4353
let
	neg_adjmatrix = Tropical.(-[Inf 5 3 2 Inf;
		5 Inf 7 Inf Inf;
		3 7 Inf 2 Inf;
		2 Inf 2 Inf 3;
		Inf Inf Inf 3 Inf])

	Compose.set_default_graphic_size(14cm, 7cm)
	sq = nodestyle(:circle, fill("white"), stroke("black"); r=0.08)
	wb = nodestyle(:circle, fill("white"); r=0.05)
	eb = bondstyle(:line, stroke("black"), linewidth(0.5mm))
	tb = textstyle(:default, fontsize(25px))
	x0 = 0.15
	x1 = 0.65
	y0 = 0.35
	y1 = 0.8
	x3 = 0.9
	y3 = 0.1
	a = (x0, y0)
	b = (x0, y1)
	c = (x1, y1)
	d = (x1, y0)
	e = (x3, y3)
	nodes = [a, b, c, d, e]
	img = canvas() do
		for (loc, label) in zip(nodes, ["A", "B", "C", "D", "E"])
			sq >> loc
			tb >> (loc, label)
		end
		for edge in [(1,2), (2,3), (3,4), (1,4), (4,5), (1,3)]
			x, y, = nodes[edge[1]],nodes[edge[2]]
			eb >> (x, y)
			wb >> ((x .+ y) ./ 2)
			tb >> ((x .+ y) ./ 2, string(-Int(neg_adjmatrix[edge...].n)))
		end
	end
	Compose.compose(context(.38, 0, .5, 1), img)
end

# ╔═╡ 25d64fd4-7b04-11eb-105a-71f98c236ec6
md"What is the shorted parth from `B` to `E`? Its adjacency matrix is given by"

# ╔═╡ 71ae3014-7b01-11eb-1707-2f6d249d8cf6
md"""
```math
A = \left[\begin{matrix}\infty & 5 & 3 & 2 & \infty\\
		5 & \infty & 7 & \infty & \infty\\
		3 & 7 & \infty & 2 & \infty\\
		2 & \infty & 2 & \infty & 3\\
		\infty & \infty & \infty & 3 & \infty\\
\end{matrix}\right]
```
"""

# ╔═╡ 0e891c90-7b0d-11eb-2616-83f1e9a21ae2
md"First, represent the negative adjacency matrix ``-A`` with tropical matrix.
"

# ╔═╡ 14b93998-7b25-11eb-22cc-f1a34c46570f
md"
The power of ``-A`` is
```math
((-A)^m)_{i_1i_m} := \max_{i_1,i_2,\ldots i_m} \left((-A)_{i_1i_2} + (-A)_{i_2i_3} +\ldots (-A)_{i_{m-1}i_m}\right)
```
Then the shortest path between ``i_1`` and ``i_m`` can be represented as
```math
-\max((-A)_{i_1i_m}, (-A)^2_{i_1i_m}, \ldots (-A)^{n-1}_{i_1i_m})
```
where ``n`` is the number of vertices. The power stops at order ``n-1`` because it shortest path contains at most ``n-1`` segments.
"

# ╔═╡ 48749c00-7b25-11eb-19f6-41b57a9808ab
md"In this example, the shortest path between `B` and `E` is 10."

# ╔═╡ 5da00572-7b04-11eb-01f1-cd94579b478e
md"
### References

* [Tropical Arithmetic and Shortest Paths](https://personalpages.manchester.ac.uk/staff/mark.muldoon/Teaching/DiscreteMaths/LectureNotes/TropicalShortestPaths.pdf)

* [Methods and Applications of (max,+) Linear Algebra](https://link.springer.com/chapter/10.1007/BFb0023465)"

# ╔═╡ 211911da-7a18-11eb-12d4-65b0dec4b8dc
md"## Tropical tensor networks
```math
\begin{align}
\cancel{Y[n] := \sum_{i,j,k,l,m} A[i,l] \times B[i,j] \times C[j,k,n] \times D[k,l,m] \times E[m]}\\

Y[n] := \max_{i,j,k,l,m} (A[i,l] + B[i,j] + C[j,k,n] + D[k,l,m] + E[m])
\end{align}
```
"

# ╔═╡ 281a5cf0-7a18-11eb-3385-c3e64f41e4da
let
	Compose.set_default_graphic_size(14cm, 7cm)
	sq = nodestyle(:circle; r=0.085)
	wb = nodestyle(:circle, fill("white"); r=0.04)
	eb = bondstyle(:line)
	tb = textstyle(:default, fontsize(25px))
	tb2 = textstyle(:default, fontsize(30px), fill("white"))
	x0 = 0.15
	x1 = 0.65
	y0 = 0.35
	y1 = 0.8
	x3 = 0.9
	y3 = 0.1
	a = (x0, y0)
	b = (x0, y1)
	c = (x1, y1)
	d = (x1, y0)
	e = (x3, y3)
	img = canvas() do
		for (loc, label) in [(a, "A"), (b, "B"), (c, "C"), (d, "D"), (e, "E")]
			sq >> loc
			tb2 >> (loc, label)
		end
		for (edge, label) in [((a, b), "i"), ((b, c), "j"), ((c, d), "k"), ((a, d), "l"), ((d,e), "m"), ((c, (0.9, 0.55)), "n")]
			eb >> edge
			wb >> ((edge[1] .+ edge[2]) ./ 2)
			tb >> ((edge[1] .+ edge[2]) ./ 2, label)
		end
	end
	Compose.compose(context(.38, 0, .5, 1), img)
end


# ╔═╡ 170515ba-7d3b-11eb-281d-457606b3f401
"black nodes are tensors"

# ╔═╡ 1e164b58-7d3b-11eb-13a0-853fb1aafaa4
"each index appear twice"

# ╔═╡ 2c47b692-7a18-11eb-2878-2976435507f5
md"""
* a tropical tensor network enumerate over internal degrees of freedom, and pick the configuration that maximizes the summation.
"""

# ╔═╡ 3237e33e-7a17-11eb-2869-b92d0801bc6e
md"
#### Reference
A Practical Introduction to Tensor Networks: Matrix Product States and Projected Entangled Pair States

[arXiv: 1306.2164](https://arxiv.org/abs/1306.2164)"

# ╔═╡ 442bcb3c-7940-11eb-18e5-d3158b74b1dc
html"""
<h1>Mapping hard problems to Tropical Tensor networks</h1>
<p>The following book covers all problems mentioned in this notebook.</p>
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
* decision: NP-Complete
* counting: #P-complete

Ising spin glass is an energy model defined on a graph

```math
-E = \sum\limits_{i,j\in E} J_{ij} s_i s_j + \sum\limits_{i\in V} h_{i} s_i
```
where $s_i,s_j,\ldots s_n \in \{-1,1\}$.
"""

# ╔═╡ d3b4f162-7ad4-11eb-271c-677cb307c447
md"## Example"

# ╔═╡ 22118a36-7a36-11eb-18c3-dd2adac6118b
md"Let's take $h_i=0$ from now on, it will not change the complexity of the problem."

# ╔═╡ b3b986aa-7a1d-11eb-17d4-e5675015b221
let
	Compose.set_default_graphic_size(14cm, 7cm)
	sq = nodestyle(:circle, fill("white"), stroke("black"); r=0.08)
	wb = nodestyle(:circle, fill("white"); r=0.05)
	eb = bondstyle(:line, stroke("red"), linewidth(0.5mm))
	eb_ = bondstyle(:line, stroke("blue"), linewidth(0.5mm))
	tb = textstyle(:default, fontsize(25px), font("times"))
	tb2 = textstyle(:default, fontsize(30px), fill("black"), font("times"))
	x0 = 0.15
	x1 = 0.65
	y0 = 0.35
	y1 = 0.8
	x3 = 0.9
	y3 = 0.1
	a = (x0, y0)
	b = (x0, y1)
	c = (x1, y1)
	d = (x1, y0)
	e = (x3, y3)
	img = canvas() do
		for (loc, label) in [(a, "s₁"), (b, "s₂"), (c, "s₃"), (d, "s₄"), (e, "s₅")]
			sq >> loc
			tb2 >> (loc, label)
		end
		for (edge, J, label) in [((a, b), 1, "J₁₂"), ((b, c), -1, "J₂₃"), ((c, d), -1, "J₃₄"), ((a, d), 1, "J₁₄"), ((d,e), 1, "J₄₅"), ((c, a), 1, "J₁₃")]
			(J > 0 ? eb : eb_) >> edge
			wb >> ((edge[1] .+ edge[2]) ./ 2)
			tb >> ((edge[1] .+ edge[2]) ./ 2, label)
		end
	end
	Compose.compose(context(.38, 0, .5, 1), img)
end

# ╔═╡ b5cd769e-7a1e-11eb-1d82-e1c265dfdd52
md"The goal is to find an optimal assignment of $s_1, s_2,\ldots s_n$ that minimizes the energy"

# ╔═╡ b52ead96-7a2a-11eb-334f-e5e5ff5867e3
md"""
```math
T_{e}(J_{ij})_{s_i s_j} = \begin{bmatrix}J_{ij} & -J_{ij} \\-J_{ij} & J_{ij}\end{bmatrix}_{s_is_j}
```
"""

# ╔═╡ 6e1c507e-7a1a-11eb-05bc-dbccc3aebdf9
md"""
```math
\begin{align}
-E = &\sum\limits_{i,j\in E} J_{ij} s_i s_j\\
=&\max_{s_1,s_2,s_3,s_4, s_5}
T_e(J_{12})_{s_1s_2}+
T_e(J_{23})_{s_2s_3}+
T_e(J_{34})_{s_3s_4}+\\
&T_e(J_{14})_{s_1s_4}+
T_e(J_{45})_{s_4s_5}+
T_e(J_{13})_{s_1s_3}
\end{align}
```
"""

# ╔═╡ 05109d30-7a29-11eb-320a-fb0b0d8e2632
md"where a spin $s_i = 1$ (or $s_i=-1$) is equivalent to $1$ (or $2$) when used in matrix indexing. 
The boond tensors are defined as"

# ╔═╡ fb5d95ea-7d25-11eb-1f17-b3936a33779a
md"Isn't this familiar? But it is not a tensor network.

```math
\begin{align}
Y[n¹] := &\max_{i,j,k,l,m} (A[i¹,l¹] + B[i²,j¹] \\
&+ C[j²,k¹,n²] + D[k²,l²,m¹] + E[m²])
\end{align}
```"

# ╔═╡ 9ad76a2e-7d26-11eb-2fab-55b656bf042a
md"A label in a tensor network appears exactly twice!"

# ╔═╡ 00ed185e-7a2d-11eb-1b27-cb834e75e916
md"## This is Einsum"

# ╔═╡ c4c87cdc-7d3b-11eb-1831-930fe8626f9e
comment("move to end")

# ╔═╡ 5efee244-7a2d-11eb-3782-b9d55086d623
md"`einsum` is a generalization of tensor networks, it allows a same indices appearing for an arbituary times. Its graphical representation is a hypergraph rather than a simple graph."

# ╔═╡ d0ecd3f2-7a2d-11eb-126d-7dab740d8e1f
let
	Compose.set_default_graphic_size(14cm, 7cm)
	sq = nodestyle(:circle, fill("black"); r=0.02)
	wb = nodestyle(:circle, fill("black"); r=0.06)
	eb = bondstyle(:line)
	tb = textstyle(:default, fontsize(20px), fill("white"))
	tb2 = textstyle(:default, fontsize(15px), fill("black"), font("times"))
	tb3 = textstyle(:default, fontsize(15px), fill("green"))
	x0 = 0.15
	x1 = 0.65
	y0 = 0.35
	y1 = 0.8
	x3 = 0.9
	y3 = 0.1
	a = (x0, y0)
	b = (x0, y1)
	c = (x1, y1)
	d = (x1, y0)
	e = (x3, y3)
	img = canvas() do
		for (edge, label) in [((a, b), "+1"), ((b, c), "-1"), ((c, d), "-1"), ((a, d), "+1"), ((d,e), "+1"), ((c, a), "+1")]
			eb >> edge
			wb >> ((edge[1] .+ edge[2]) ./ 2)
			tb >> ((edge[1] .+ edge[2]) ./ 2, label == "+1" ? "Tₑ" : "Tₑ'")
		end
		for (x, deg) in zip([a .- (0,0.05), b .+ (0,0.05), c.+ (0,0.05), d .+ (0.15, 0), e .- (0,0.05)], [3,2,3,3,1])
			tb2 >> (x, "degree $deg")
		end
	end
	Compose.compose(context(.38, 0, .5, 1), img)
end

# ╔═╡ 8d8865c0-7d27-11eb-0380-8f0a6af81e92
md"In the above graph, there are three hyperedges of size 3, one hyperedge of size 2 and one hyperedge of size 1."

# ╔═╡ 5bcae146-7d27-11eb-2ad7-a7e419ec49ae
md"## Converting Einsum to Tensor networks
einsum can be converted to tensor network by introducing δ tensors"

# ╔═╡ 96290770-7a20-11eb-0ac8-33a6492c7b12
md"To map the spin glass problem to a tropical tensor network, we place a tensor $T_v^d$ of rank-d at each vertex, where $d$ is the degree of the vertex."

# ╔═╡ c1f90d6c-7a1d-11eb-2843-f971b5f6f3b0
md"""
```math
T_{v}^{n} = \begin{cases}
 h, & s_i = s_j =\ldots s_n = 0\\
 -h, & s_i = s_j =\ldots s_n = 1\\
 -\infty, &otherwise
\end{cases}
```
"""

# ╔═╡ 64e08a56-7a36-11eb-29fd-03662b4d6612
md"for $h=0$, this is equivalent to the δ tensor."

# ╔═╡ 7753f426-7d28-11eb-01af-0599499339ae
md"""
```math
\begin{align}
-E = &\max_{s_1,s_2,s_3,s_4, s_5}
T_e(J_{12})_{s_1s_2}+
T_e(J_{23})_{s_2's_3}+
T_e(J_{34})_{s_3s_4}+\\
&T_e(J_{14})_{s_1's_4'}+
T_e(J_{45})_{s_4''s_5}+
T_e(J_{13})_{s_1''s_3'}+\\
&(T_v^2)_{s_1s_1's_1''}+
(T_v^3)_{s_2s_2'}+
(T_v^3)_{s_3s_3's_3''}+\\
&(T_v^1)_{s_4s_4's_4''}+
(T_v^1)_{s_5}
\end{align}
```
"""

# ╔═╡ f54119ca-7a1e-11eb-1bec-bf855e34658d
let
	Compose.set_default_graphic_size(14cm, 7cm)
	sq = nodestyle(:circle, fill("black"); r=0.02)
	wb = nodestyle(:circle, fill("black"); r=0.06)
	eb = bondstyle(:line)
	tb = textstyle(:default, fontsize(20px), fill("white"))
	tb2 = textstyle(:default, fontsize(20px), fill("black"), font("times"))
	tb3 = textstyle(:default, fontsize(15px), fill("green"))
	x0 = 0.15
	x1 = 0.65
	y0 = 0.35
	y1 = 0.8
	x3 = 0.9
	y3 = 0.1
	a = (x0, y0)
	b = (x0, y1)
	c = (x1, y1)
	d = (x1, y0)
	e = (x3, y3)
	img = canvas() do
		for (loc, label) in [(a, "↑"), (b, "↓"), (c, "↑"), (d, "↑"), (e, "↓")]
			sq >> loc
			#tb2 >> (loc, label)
		end
		for (edge, label) in [((a, b), "+1"), ((b, c), "-1"), ((c, d), "-1"), ((a, d), "+1"), ((d,e), "+1"), ((c, a), "+1")]
			eb >> edge
			wb >> ((edge[1] .+ edge[2]) ./ 2)
			tb >> ((edge[1] .+ edge[2]) ./ 2, label == "+1" ? "Tₑ" : "Tₑ'")
		end
		for lt in [((0.08, 0.3),"Tᵥ³"), ((0.08,0.75), "Tᵥ²"), ((0.75, 0.35),"Tᵥ³"), ((0.75,0.8), "Tᵥ²"), ((0.97,0.1), "Tᵥ¹")]
			tb2 >> lt
		end
		for (i,edge) in enumerate([(a, b), (b, c), (c, d), (a, d), (d,e), (c, a)])
			#tb3 >> (edge[1] .* 0.2 .+ edge[2] .* 0.8, "($i, T)")
			#tb3 >> (edge[1] .* 0.8 .+ edge[2] .* 0.2, "($i, F)")
		end
	end
	Compose.compose(context(.38, 0, .5, 1), img)
end

# ╔═╡ 37102544-7abf-11eb-3ac4-6702dfc55425
html"""big circles are for vertex tensors and edge tensors respectively.</p>"""

# ╔═╡ 695e405c-786d-11eb-0a6e-bb776d9626ad
md"
## Using Tropical numbers for counting
We define the multiplication and addition of the tuple: $(x_1, n_1) \odot (x_2,n_2) = (x_1 + x_2, n_1\cdot n_2)$ and 
```math
\begin{equation}
    (x_1, n_1)\oplus (x_2, n_2) = \begin{cases}
 (x_1\oplus x_2, \, n_1 + n_2 ) & \text{if $x_1 = x_2$} \\
 (x_1\oplus x_2,\, n_1 ) & \text{if $x_1>x_2$} \\
 (x_1\oplus x_2,\, n_2 )& \text{if $x_1 < x_2$}
 \end{cases}.
\end{equation}
```

Where the ``x`` is the value and ``n`` is the counting.
zero elements are ``(-\infty, 1)`` and one elements are ``(0, 0)``
"

# ╔═╡ 875e8632-7d3d-11eb-1472-f9d949d14a45
comment("why is this? image->two graphs")

# ╔═╡ 4190393a-7ac4-11eb-3ac6-eb8e3574fdc9
md"There are 2 degenerate ground states with energy -4, they are"

# ╔═╡ 56fdb22c-7ac4-11eb-2831-a777d9ca89f3
let
	Compose.set_default_graphic_size(16cm, 7cm)
	sq = nodestyle(:circle, fill("white"), stroke("black"); r=0.08)
	eb = bondstyle(:line, stroke("red"), linewidth(0.5mm))
	eb_ = bondstyle(:line, stroke("blue"), linewidth(0.5mm))
	tb = textstyle(:default, fontsize(25px), font("times"))
	tb2 = textstyle(:default, fontsize(30px), fill("black"), font("times"))
	x0 = 0.15
	x1 = 0.65
	y0 = 0.35
	y1 = 0.8
	x3 = 0.9
	y3 = 0.1
	a = (x0, y0)
	b = (x0, y1)
	c = (x1, y1)
	d = (x1, y0)
	e = (x3, y3)
	img1 = canvas() do
		for (loc, label) in [(a, "↓"), (b, "↓"), (c, "↑"), (d, "↓"), (e, "↓")]
			sq >> loc
			tb2 >> (loc, label)
		end
		for (edge, J, label) in [((a, b), 1, "J₁₂"), ((b, c), -1, "J₂₃"), ((c, d), -1, "J₃₄"), ((a, d), 1, "J₁₄"), ((d,e), 1, "J₄₅"), ((c, a), 1, "J₁₃")]
			(J > 0 ? eb : eb_) >> edge
		end
	end
	img2 = canvas() do
		for (loc, label) in [(a, "↑"), (b, "↑"), (c, "↓"), (d, "↑"), (e, "↑")]
			sq >> loc
			tb2 >> (loc, label)
		end
		for (edge, J, label) in [((a, b), 1, "J₁₂"), ((b, c), -1, "J₂₃"), ((c, d), -1, "J₃₄"), ((a, d), 1, "J₁₄"), ((d,e), 1, "J₄₅"), ((c, a), 1, "J₁₃")]
			(J > 0 ? eb : eb_) >> edge
		end
	end
	Compose.compose(context(), (context(.05, 0, .5*7/8, 1), img1),
		(context(.55, 0, .5*7/8, 1), img2))
end

# ╔═╡ e739e74c-7af0-11eb-104f-5f94da1bf0be
md"## In our paper"

# ╔═╡ 80d6c2b6-7aef-11eb-1bf5-5d4f266dfa73
let
	Compose.set_default_graphic_size(7cm, 7cm)
	img = viz(SquareLattice(18,18); node_style=nodestyle(:circle; r=0.012), text_style=textstyle(:default, fill("transparent")))
	leftright(updown(
		img,
		md"Ground-state energy, entropy, and computational
time of ``q = 3`` state Potts spin glass model on square lattices of sizes ``n = 4\sim 18``. Each data point is
averaged over ``100`` random instances computed on a single
GPU. As a comparison, the existing branch-and-cut method
with the Semi-Definition Programming energy lower bounds
method on the same model works up to ``9 \times 9`` lattices
(using 10 hours)"
		; width=300)
		,html"""<img src="https://user-images.githubusercontent.com/6257240/109578583-6d8c7680-7ac5-11eb-93eb-1b2748f2c90b.png" width=250px/>""")
end

# ╔═╡ 344042b4-793d-11eb-3d6f-43eb2a4db9f4
md"""
## Maximum independent set
* hardness: NP-Complete
* counting: #P-complete
"""

# ╔═╡ 80d764a8-7afd-11eb-3fb8-79169ca56c7e
md"A maximal independent set (MIS) is an independent set that is not a subset of any other independent set. Let's still use this graph as an example."

# ╔═╡ 04d11828-7afa-11eb-3e73-1bbecf566f74
let
	Compose.set_default_graphic_size(14cm, 7cm)
	sq = nodestyle(:circle, fill("white"), stroke("black"); r=0.08)
	wb = nodestyle(:circle, fill("white"); r=0.05)
	eb = bondstyle(:line, stroke("black"), linewidth(0.5mm))
	tb = textstyle(:default, fontsize(25px), font("times"))
	tb2 = textstyle(:default, fontsize(30px), fill("black"), font("times"))
	x0 = 0.15
	x1 = 0.65
	y0 = 0.35
	y1 = 0.8
	x3 = 0.9
	y3 = 0.1
	a = (x0, y0)
	b = (x0, y1)
	c = (x1, y1)
	d = (x1, y0)
	e = (x3, y3)
	img = canvas() do
		for (loc, label) in [(a, "s₁"), (b, "s₂"), (c, "s₃"), (d, "s₄"), (e, "s₅")]
			sq >> loc
			tb2 >> (loc, label)
		end
		for (edge, J, label) in [((a, b), 1, "J₁₂"), ((b, c), -1, "J₂₃"), ((c, d), -1, "J₃₄"), ((a, d), 1, "J₁₄"), ((d,e), 1, "J₄₅"), ((c, a), 1, "J₁₃")]
			eb >> edge
			#wb >> ((edge[1] .+ edge[2]) ./ 2)
			#tb >> ((edge[1] .+ edge[2]) ./ 2, label)
		end
	end
	Compose.compose(context(.38, 0, .5, 1), img)
end

# ╔═╡ 5d95a598-7afa-11eb-10eb-db79fa44dd2a
md"""
There are different ways to map a MIS problem to tensor networks. 
"""

# ╔═╡ 14617e4e-7d37-11eb-05d0-0f223b0f8721
md"## Mapping No.1"

# ╔═╡ d29470d4-7afa-11eb-0afc-a34e39d49aa5
let
	Compose.set_default_graphic_size(14cm, 7cm)
	sq = nodestyle(:circle, fill("black"); r=0.02)
	wb = nodestyle(:circle, fill("black"); r=0.06)
	eb = bondstyle(:line)
	tb = textstyle(:default, fontsize(20px), fill("white"))
	tb2 = textstyle(:default, fontsize(20px), fill("black"), font("times"))
	tb3 = textstyle(:default, fontsize(15px), fill("green"))
	x0 = 0.15
	x1 = 0.65
	y0 = 0.35
	y1 = 0.8
	x3 = 0.9
	y3 = 0.1
	a = (x0, y0)
	b = (x0, y1)
	c = (x1, y1)
	d = (x1, y0)
	e = (x3, y3)
	img = canvas() do
		for (loc, label) in [(a, "↑"), (b, "↓"), (c, "↑"), (d, "↑"), (e, "↓")]
			sq >> loc
			#tb2 >> (loc, label)
		end
		for (edge, label) in [((a, b), "+1"), ((b, c), "-1"), ((c, d), "-1"), ((a, d), "+1"), ((d,e), "+1"), ((c, a), "+1")]
			eb >> edge
			wb >> ((edge[1] .+ edge[2]) ./ 2)
			tb >> ((edge[1] .+ edge[2]) ./ 2, label == "+1" ? "Tₑ" : "Tₑ'")
		end
		for lt in [((0.08, 0.3),"Tᵥ³"), ((0.08,0.75), "Tᵥ²"), ((0.75, 0.35),"Tᵥ³"), ((0.75,0.8), "Tᵥ²"), ((0.97,0.1), "Tᵥ¹")]
			tb2 >> lt
		end
		for (i,edge) in enumerate([(a, b), (b, c), (c, d), (a, d), (d,e), (c, a)])
			tb3 >> (edge[1] .* 0.2 .+ edge[2] .* 0.8, "($i, T)")
			tb3 >> (edge[1] .* 0.8 .+ edge[2] .* 0.2, "($i, F)")
		end
	end
	Compose.compose(context(.38, 0, .5, 1), img)
end

# ╔═╡ 75c37046-7b1b-11eb-00f5-7fc49f73f4d9
md"1. We can put a bond matrix at the bond to describe the independence restriction
```math
T_{e} = \begin{bmatrix}0 & 0 \\0 & -\infty\end{bmatrix}
```

The vertex tensors are for counting the number of vertices
```math
(T_{v}^{n})_{s_i s_j} = \begin{cases}
 0, & s_i = s_j =\ldots s_n = 0\\
 1, & s_i = s_j =\ldots s_n = 1\\
 -\infty, &otherwise
\end{cases}
```
where $s_i,s_j,\ldots s_n \in \{0,1\}$."

# ╔═╡ 28534a16-7d37-11eb-0747-316a35af4f73
md"## Mapping No.2"

# ╔═╡ c4fc59c6-7d3e-11eb-1267-239a3240519c
comment("move to the end")

# ╔═╡ 891c39b8-7b1b-11eb-2d70-67fd24021027
md"2. Equivalently, we can transfer the weights on the vertex to bond tensors so that the vertex tensors can become δ tensors.
```math
T_{e} = \begin{bmatrix}0 & \frac{1}{w_{j→e}} \\\frac{1}{w_{i→e}} & -\infty\end{bmatrix}
```

where $i,j$ are source and destination node of $e$, $w_{i→e}$ is the weight transfered from vertex $i$ to $e$ that $\sum\limits_e w_{i\rightarrow e}=1$. Then vertex tensor now is
```math
T_{v}^{n} = \delta^{n,q=2}
```
"

# ╔═╡ 3c06092c-7d37-11eb-01f5-77b6cb8967eb
md"## Some Results
* a comparison with MoMC
* sparse tensor contraction for larger size
* graph instances (high degeneraacy - low degeneracy)
"

# ╔═╡ 5caa1b7c-7b21-11eb-394d-379351fe5170
md"## A new solution to the SSSS challenge"

# ╔═╡ 1c4b19d2-7b30-11eb-007b-ab03052b22d2
md"With the greedy contraction order"

# ╔═╡ f14f7fa8-7d31-11eb-2f4d-2bc1ddbacdb6
md"number of steps = $(@bind nstep_c60 Slider(0:150; show_value=true, default=29))"

# ╔═╡ 4c137484-7b30-11eb-2fb1-190d8beebbc3
md"Finding the optimal contraction order is know as NP-hard, however, there are [some heuristic algorithms](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.060503) can find good contraction path in limited time.
"

# ╔═╡ e302bd1c-7ab5-11eb-03f6-69dcbb817354
md"# Future Directions

* Finding the optimal contracting order of `einsum`,
    - einsum has less redundancy, with less effort in mapping, and potentially faster in contraction,
    - tensor network can utilize BLAS easily,
    - the contraction order for tensor networks are more studies.
    - the optimal configuration
* TropicalBLAS,
    - TropicalGEMM on CPU for Tropical numbers (DONE),
    - TropicalGEMM on GPU (WIP),
    - TropicalGEMM on CPU for CountingTropical numbers,
"

# ╔═╡ d53a5e5a-7b33-11eb-0a65-4d4844e9cf0a
md"# Utilities"

# ╔═╡ e3e88ea6-7d1b-11eb-2400-17565c85cac2
begin
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
	
	function mis_bondtensor(::Type{T}) where T
		res = ones(T, 2, 2)
		res[2, 2] = zero(T)
		return res
	end
	
	function mis_vertextensor(::Type{T}, n::Int) where T
		res = zeros(T, fill(2, n)...)
		res[1] = one(T)
		res[end] = T(1)
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

	function δtensor(::Type{T}, q::Int, n::Int) where {T}
		res = zeros(T, fill(q, n)...)
		for i=1:q
			res[fill(i, n)...] = one(T)
		end
		res
	end
	potts_vertextensor(args...) = δtensor(args...)
	
	function build_tensornetwork(; vertices, vertex_arrays, edges, edge_arrays)
		TensorNetwork([
		# vertex tensors
		[LabeledTensor(vertex_arrays[i], [(j, v==e[1]) for (j, e) in enumerate(edges) if v ∈ e]) for (i, v) in enumerate(vertices)]...,
		# bond tensors
		[LabeledTensor(edge_arrays[j], [(j, true), (j, false)]) for j=1:length(edges)]...
	])
	end
end;

# ╔═╡ 9b1dc21a-7896-11eb-21f6-bfe9b4dc9ccf
let
	θ2 = 0.5
	ϕ2 = 0.8
	Compose.set_default_graphic_size(12cm, 10cm)
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
	Compose.compose(context(),
		(context(), text(0.75, 0.75, "60 vertices")),
		(context(), text(0.75, 0.8, "120 edges")),
		(context(0.4,0.5, 1.0, 1.25), fig),
		)
end

# ╔═╡ c1c74e70-7b2c-11eb-2f26-21f54ad00fb2
let
	function contraction_mask(tnet, tree)
		contraction_mask!(tnet, tree, [zeros(Bool, length(tnet))])
	end
	function contraction_mask!(tnet, tree, results)
		if tree isa Integer
			res = copy(results[end])
			@assert res[tree] == false
			res[tree] = true
			push!(results, res)
		else
			contraction_mask!(tnet, tree.left, results)
			contraction_mask!(tnet, tree.right, results)
		end
		return results
	end

	c60_xy = fullerene();
	c60_edges = [(i=>j) for (i,(i2,j2,k2)) in enumerate(c60_xy), (j,(i1,j1,k1)) in enumerate(c60_xy) if i<j && (i2-i1)^2+(j2-j1)^2+(k2-k1)^2 < 5.0];
	
	c60_tnet = let
		T = CountingTropical{Float64}
		build_tensornetwork(
			vertices=1:60,
			vertex_arrays = [ising_vertextensor(T, 3, 0.0) for j=1:length(c60_xy)],
			edges = c60_edges,
			edge_arrays = [ising_bondtensor(T, -1.0) for i = 1:length(c60_edges)]
		)
	end;
	tcs, scs, c60_trees = (Random.seed!(2); trees_greedy(c60_tnet; strategy="min_reduce"));
	c60_contraction_masks = contraction_mask(c60_tnet, c60_trees[])

	#SimpleTensorNetworks.contract(c60_tnet, c60_trees[]).array[]
	mask = c60_contraction_masks[nstep_c60+1]
	Compose.set_default_graphic_size(12cm, 12cm)
	cam_position = SVector(0.0, 0.0, 0.5)
	rot = RotY(0.5)*RotX(0.8)
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
	Compose.compose(context(0.5,0.35, 1.0, 1.0), fig)
end

# ╔═╡ Cell order:
# ╟─c456b902-7959-11eb-03ba-dd14a2cd5758
# ╟─7fa86de2-7ac5-11eb-3d55-c9b6cb57c5d6
# ╟─5bb40ad6-7b33-11eb-0b31-63d5e47fa0e7
# ╟─121b4926-7aba-11eb-30e1-7b8edd4f0166
# ╟─265649e4-7abb-11eb-1a02-a3c101cef89d
# ╟─7bdf517e-79ff-11eb-38a3-49c02d94d943
# ╟─9b1dc21a-7896-11eb-21f6-bfe9b4dc9ccf
# ╟─88e14ef2-7af1-11eb-23d6-b34b1eff8f87
# ╟─a387b018-79ff-11eb-0383-1fcf82853afc
# ╟─3205a536-7a17-11eb-3473-b71305c96ca4
# ╟─be76e52a-7852-11eb-179b-afbc6efcab55
# ╟─d0b54b76-7852-11eb-2398-0911380fa090
# ╟─9dbbcb08-7d40-11eb-3d58-d1ccf9a9d803
# ╟─90ee0580-7d38-11eb-36d3-299cefec4ddf
# ╠═908952d0-7d0d-11eb-2cf7-6fa6d1fbee88
# ╟─98ae0960-797d-11eb-3646-c5b7e05d3f7c
# ╟─86921d00-7a17-11eb-2695-add5f9eeda5b
# ╟─915e8096-7a17-11eb-177d-a39ffed7ca91
# ╠═1cbb86e4-7d39-11eb-05d9-c3ab43a03298
# ╟─958c489e-7a17-11eb-2be4-cfdb44da7d2f
# ╟─ddca618c-7d22-11eb-0c19-e1a149e19f3c
# ╟─3d7ca08c-7b01-11eb-1d78-af35dc7e577c
# ╠═4fbec178-7d39-11eb-1f27-13ac45f53e45
# ╟─742fdf74-7d0e-11eb-026f-81ee2c7c4353
# ╟─25d64fd4-7b04-11eb-105a-71f98c236ec6
# ╟─71ae3014-7b01-11eb-1707-2f6d249d8cf6
# ╟─0e891c90-7b0d-11eb-2616-83f1e9a21ae2
# ╟─14b93998-7b25-11eb-22cc-f1a34c46570f
# ╟─48749c00-7b25-11eb-19f6-41b57a9808ab
# ╟─5da00572-7b04-11eb-01f1-cd94579b478e
# ╟─211911da-7a18-11eb-12d4-65b0dec4b8dc
# ╟─281a5cf0-7a18-11eb-3385-c3e64f41e4da
# ╠═170515ba-7d3b-11eb-281d-457606b3f401
# ╠═1e164b58-7d3b-11eb-13a0-853fb1aafaa4
# ╟─2c47b692-7a18-11eb-2878-2976435507f5
# ╟─3237e33e-7a17-11eb-2869-b92d0801bc6e
# ╟─442bcb3c-7940-11eb-18e5-d3158b74b1dc
# ╟─f7208b6e-793c-11eb-0dfa-0d63752ba53e
# ╟─d3b4f162-7ad4-11eb-271c-677cb307c447
# ╟─22118a36-7a36-11eb-18c3-dd2adac6118b
# ╟─b3b986aa-7a1d-11eb-17d4-e5675015b221
# ╟─b5cd769e-7a1e-11eb-1d82-e1c265dfdd52
# ╟─b52ead96-7a2a-11eb-334f-e5e5ff5867e3
# ╟─6e1c507e-7a1a-11eb-05bc-dbccc3aebdf9
# ╟─05109d30-7a29-11eb-320a-fb0b0d8e2632
# ╟─fb5d95ea-7d25-11eb-1f17-b3936a33779a
# ╟─9ad76a2e-7d26-11eb-2fab-55b656bf042a
# ╟─00ed185e-7a2d-11eb-1b27-cb834e75e916
# ╠═c4c87cdc-7d3b-11eb-1831-930fe8626f9e
# ╟─5efee244-7a2d-11eb-3782-b9d55086d623
# ╟─d0ecd3f2-7a2d-11eb-126d-7dab740d8e1f
# ╟─8d8865c0-7d27-11eb-0380-8f0a6af81e92
# ╟─5bcae146-7d27-11eb-2ad7-a7e419ec49ae
# ╟─96290770-7a20-11eb-0ac8-33a6492c7b12
# ╟─c1f90d6c-7a1d-11eb-2843-f971b5f6f3b0
# ╟─64e08a56-7a36-11eb-29fd-03662b4d6612
# ╟─7753f426-7d28-11eb-01af-0599499339ae
# ╟─f54119ca-7a1e-11eb-1bec-bf855e34658d
# ╟─37102544-7abf-11eb-3ac4-6702dfc55425
# ╟─695e405c-786d-11eb-0a6e-bb776d9626ad
# ╠═875e8632-7d3d-11eb-1472-f9d949d14a45
# ╟─4190393a-7ac4-11eb-3ac6-eb8e3574fdc9
# ╟─56fdb22c-7ac4-11eb-2831-a777d9ca89f3
# ╟─e739e74c-7af0-11eb-104f-5f94da1bf0be
# ╟─80d6c2b6-7aef-11eb-1bf5-5d4f266dfa73
# ╟─344042b4-793d-11eb-3d6f-43eb2a4db9f4
# ╟─80d764a8-7afd-11eb-3fb8-79169ca56c7e
# ╟─04d11828-7afa-11eb-3e73-1bbecf566f74
# ╟─5d95a598-7afa-11eb-10eb-db79fa44dd2a
# ╟─14617e4e-7d37-11eb-05d0-0f223b0f8721
# ╟─d29470d4-7afa-11eb-0afc-a34e39d49aa5
# ╟─75c37046-7b1b-11eb-00f5-7fc49f73f4d9
# ╟─28534a16-7d37-11eb-0747-316a35af4f73
# ╠═c4fc59c6-7d3e-11eb-1267-239a3240519c
# ╟─891c39b8-7b1b-11eb-2d70-67fd24021027
# ╟─3c06092c-7d37-11eb-01f5-77b6cb8967eb
# ╟─5caa1b7c-7b21-11eb-394d-379351fe5170
# ╟─1c4b19d2-7b30-11eb-007b-ab03052b22d2
# ╟─f14f7fa8-7d31-11eb-2f4d-2bc1ddbacdb6
# ╟─c1c74e70-7b2c-11eb-2f26-21f54ad00fb2
# ╟─4c137484-7b30-11eb-2fb1-190d8beebbc3
# ╟─e302bd1c-7ab5-11eb-03f6-69dcbb817354
# ╟─d53a5e5a-7b33-11eb-0a65-4d4844e9cf0a
# ╟─e3e88ea6-7d1b-11eb-2400-17565c85cac2
