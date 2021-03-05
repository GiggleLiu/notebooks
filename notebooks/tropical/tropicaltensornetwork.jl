### A Pluto.jl notebook ###
# v0.12.18

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
	PlutoUI.TableOfContents()
end

# ╔═╡ 5bb40ad6-7b33-11eb-0b31-63d5e47fa0e7
using TropicalNumbers,  # tropical number type
		LightGraphs,	# graph operations
		Random,
    	SimpleTensorNetworks  # tensor network contraction

# ╔═╡ 1749c0f2-7a2a-11eb-1932-07a7f920b0da
using OMEinsum

# ╔═╡ dfa8834c-e8c6-49b4-8bde-0816b573cbee
html"""
<style>
body {
counter-reset: section subsection example}

h2::before {
counter-reset: subsection;
  counter-increment: section;
  content: "Sec. " counter(section) ": ";
}
"""

# ╔═╡ 121b4926-7aba-11eb-30e1-7b8edd4f0166
md"""# Tropical tensor networks for solving combinatoric optimization problems


$(HTML("<br><p><big><strong>Tropical tensor network for ground states of spin glasses</strong></big></p>
<p>Phys. Rev. Lett. (26 January 2021)</p>
<p>Jin-Guo Liu, Lei Wang, and Pan Zhang</p>"))


[arxiv 2008.06888](https://arxiv.org/abs/2008.06888)
"""

# ╔═╡ 3205a536-7a17-11eb-3473-b71305c96ca4
md"## A minimum introduction to tensor networks"

# ╔═╡ 3208fd8a-7a17-11eb-35ce-4d6b141c1aff
md"####  A graphical representation of matrix multiplication
```math
Y[i,j] := \sum_k A[i,k] \times B[k,j]
```
"

# ╔═╡ 32116a92-7a17-11eb-228f-0713510d0348
let
	Compose.set_default_graphic_size(15cm, 10/3*cm)
	sq = nodestyle(:square; r=0.08)
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


# ╔═╡ 3217d1ca-7a17-11eb-14eb-a77ccfd983a0
md"
* a matrix is a tensor of rank-2, it is represented as a vertex with two labeled legs,
* the **contraction** (a generalization of matrix multiplication to tensors) is represented by connecting legs with same labels, where connecting a pair of legs means **sum-product** over a specific degree of freedom corresponding to the leg label,
"

# ╔═╡ 3221a326-7a17-11eb-0fe6-f75798a411b9
md"#### A graphical representation of tensor networks
```math
Y[n] := \sum_{i,j,k,l,m} A[i,l] \times B[i,j] \times C[j,k,n] \times D[k,l,m] \times E[m]
```
"

# ╔═╡ 32277c3a-7a17-11eb-3763-af68dbb81465
let
	Compose.set_default_graphic_size(14cm, 7cm)
	sq = nodestyle(:square; r=0.07)
	wb = nodestyle(:square, fill("white"); r=0.04)
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

# ╔═╡ 322d2958-7a17-11eb-2deb-613b7680a5bb
md"""
* A tensor of rank ``n`` is represented as a vertex with ``n`` labeled legs,
* A tensor network is represented as an **open simple graph**. It is a representation of the contraction of multiple tensors, the output is a tensor of rank-(the number of unpaired legs)
"""

# ╔═╡ 3237e33e-7a17-11eb-2869-b92d0801bc6e
md"
#### Reference
A Practical Introduction to Tensor Networks: Matrix Product States and Projected Entangled Pair States

[arXiv: 1306.2164](https://arxiv.org/abs/1306.2164)"

# ╔═╡ ec841be8-7a16-11eb-3337-376e26b7da25
md"## Tropical numbers and Tropical Tensor networks"

# ╔═╡ be76e52a-7852-11eb-179b-afbc6efcab55
md"#### Tropical algebra"

# ╔═╡ d0b54b76-7852-11eb-2398-0911380fa090
md"""
Tropical algebra is defined by replacing the usual sum and product operators for ordinary real numbers with the max and sum operators respectively 

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

# ╔═╡ 5d16a472-785b-11eb-1b94-dd6d8f860c24
md"""
#### It is a semi-ring!

The $\oplus$ and $\odot$ operators still have commutative, associative, and distributive properties. However, since there is no additive inverse, the $\oplus$ and $\odot$ and operations define a semiring over ``\mathbb R \cup \{-\infty\}``. 
"""

# ╔═╡ 7da42d08-7ad2-11eb-095a-87aedba91b35
md"## Zeros and ones"

# ╔═╡ 3372871c-785b-11eb-3092-4bbc419cb788
md"One sees that $-\infty$ acts as zero element for the tropical number since  $-\infty \oplus x = x  $ and $-\infty \odot x = -\infty$. "

# ╔═╡ 2173a6cc-785b-11eb-1ab6-7fb875224dd9
# Float64 is the storage type of Tropical number
zero(Tropical{Float64})

# ╔═╡ 518b7d4e-785b-11eb-3b7c-1389065b9cbd
md"
On the other hand, $0$ acts as the multiplicative identity since $0 \odot x = x$"

# ╔═╡ 2868b292-785b-11eb-015e-6b5613bd9e39
one(Tropical{Float64})

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


# ╔═╡ 958c489e-7a17-11eb-2be4-cfdb44da7d2f
md"
* to avoid confusion, we use **circles** to represent tropical tensors in this notebook,
* connecting a pair of legs means **max-sum** over a specific degree of freedom corresponding to the label,
"

# ╔═╡ 3d7ca08c-7b01-11eb-1d78-af35dc7e577c
md"## Example: Tropical matrix multiplication to find the shortest path"

# ╔═╡ 25d64fd4-7b04-11eb-105a-71f98c236ec6
md"What is the shorted parth from `B` to `E`?"

# ╔═╡ 0e891c90-7b0d-11eb-2616-83f1e9a21ae2
md"First, represent the negative adjacency matrix ``-A`` with tropical matrix.
"

# ╔═╡ 71ae3014-7b01-11eb-1707-2f6d249d8cf6
# we broadcast the `Tropical` function to each element of the adjacency matrix
neg_adjmatrix = Tropical.(-[Inf 5 3 2 Inf;
		5 Inf 7 Inf Inf;
		3 7 Inf 2 Inf;
		2 Inf 2 Inf 3;
		Inf Inf Inf 3 Inf])

# ╔═╡ 368c9728-7b01-11eb-3c30-c9d8b4d00ace
let
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

# ╔═╡ 14b93998-7b25-11eb-22cc-f1a34c46570f
md"
The power of ``-A`` is
```math
((-A)^m)_{i_1i_m} := \max_{i_1,i_2,\ldots i_m} \left((-A)_{i_1i_2} + (-A)_{i_2i_3} +\ldots (-A)_{i_{m-1}i_m}\right)
```
Then the shortest path between ``i_1`` and ``i_m`` can be represented as
```math
-\max((-A)_{i_1i_m}, (-A)^2_{i_1i_m}, \ldots (-A)^n_{i_1i_m})
```
where ``n`` is the number of vertices. The power stops at order ``n`` because it shortest path contains at most ``n`` edges.
"

# ╔═╡ 48749c00-7b25-11eb-19f6-41b57a9808ab
md"Hence, the shortest path between `B` and `E` is"

# ╔═╡ 8bd6eee8-7b03-11eb-2a79-694e7bc29bd6
# the `Tropical{Float64}` has only one field `n` of type `Float64`
# `2` and `5` are for indexing `B` and `E`
-maximum([(neg_adjmatrix^n)[2,5].n for n=1:5])

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


# ╔═╡ 2c47b692-7a18-11eb-2878-2976435507f5
md"""
* a tropical tensor network enumerate over internal degrees of freedom, and pick the configuration that maximizes the summation.
"""

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

# ╔═╡ c85217b6-7a23-11eb-04c5-fb4dc9c02ef1
J12, J23, J34, J14, J45, J13 = [1, -1, -1, 1, 1, 1]

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
		for (edge, J, label) in [((a, b), J12, "J₁₂"), ((b, c), J23, "J₂₃"), ((c, d), J34, "J₃₄"), ((a, d), J14, "J₁₄"), ((d,e), J45, "J₄₅"), ((c, a), J13, "J₁₃")]
			(J > 0 ? eb : eb_) >> edge
			wb >> ((edge[1] .+ edge[2]) ./ 2)
			tb >> ((edge[1] .+ edge[2]) ./ 2, label)
		end
	end
	Compose.compose(context(.38, 0, .5, 1), img)
end

# ╔═╡ b5cd769e-7a1e-11eb-1d82-e1c265dfdd52
md"The goal is to find an optimal assignment of $s_1, s_2,\ldots s_n$ that minimizes the energy."

# ╔═╡ 00ed185e-7a2d-11eb-1b27-cb834e75e916
md"### Mapping the problem to Einsum is straight-forward"

# ╔═╡ 9e6fbf32-7a2f-11eb-17cb-9167d6a34281
md"We can rewrite the energy definition in tropical contraction format"

# ╔═╡ 6e1c507e-7a1a-11eb-05bc-dbccc3aebdf9
md"""
```math
\begin{align}
-E = &\max_{s_1,s_2,s_3,s_4, s_5}
T_e(J_{12})_{s_1s_2}+
T_e(J_{23})_{s_2s_3}+
T_e(J_{34})_{s_3s_4}+\\
&T_e(J_{14})_{s_1s_4}+
T_e(J_{45})_{s_4s_5}+
T_e(J_{13})_{s_1s_3}
\end{align}
```
"""

# ╔═╡ b52ead96-7a2a-11eb-334f-e5e5ff5867e3
md"""
where edge tensors are defined as

```math
T_{e}(J_{ij})_{s_i s_j} = \begin{bmatrix}J_{ij} & -J_{ij} \\-J_{ij} & J_{ij}\end{bmatrix}_{s_is_j}
```
"""

# ╔═╡ 05109d30-7a29-11eb-320a-fb0b0d8e2632
md"where a spin $s_i = 1$ (or $s_i=-1$) is equivalent to $1$ (or $2$) when used in matrix indexing."

# ╔═╡ d0ecd3f2-7a2d-11eb-126d-7dab740d8e1f
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
		for (edge, label) in [((a, b), "+1"), ((b, c), "-1"), ((c, d), "-1"), ((a, d), "+1"), ((d,e), "+1"), ((c, a), "+1")]
			eb >> edge
			wb >> ((edge[1] .+ edge[2]) ./ 2)
			tb >> ((edge[1] .+ edge[2]) ./ 2, label == "+1" ? "Tₑ" : "Tₑ'")
		end
	end
	Compose.compose(context(.38, 0, .5, 1), img)
end

# ╔═╡ 9c860e2a-7a2e-11eb-231f-63e9aca1daa0
md"""This is not a tensor network since indices $s_{1}, s_3$ and $s_{4}$ appears more than twice and inner degree of freedom $s_5$ appears only once, while all indices appears once for open legs or and twice for inner degree of freedoms in a tensor network. However, we can still use `einsum` to contract this graph.
"""

# ╔═╡ 5efee244-7a2d-11eb-3782-b9d55086d623
md"`einsum` is a generalization of tensor networks, it allows a same indices appearing for an arbituary times. Its graphical representation is a hypergraph rather than a simple graph. In the above graph, there are three hyperedges of size 3, one hyperedge of size 2 and one hyperedge of size 1."

# ╔═╡ 023ebf7c-7b36-11eb-1c9f-430773395534
md"Here, the Einstein summation notation is consistent with [numpy's einsum notation](https://ajcr.net/Basic-guide-to-einsum/)"

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
			tb3 >> (edge[1] .* 0.2 .+ edge[2] .* 0.8, "($i, T)")
			tb3 >> (edge[1] .* 0.8 .+ edge[2] .* 0.2, "($i, F)")
		end
	end
	Compose.compose(context(.38, 0, .5, 1), img)
end

# ╔═╡ 37102544-7abf-11eb-3ac4-6702dfc55425
html"""<p>where the <span style="color:green">green</span> texts are labels for tensor legs. Small circles and big circles are for vertex tensors and edge tensors respectively.</p>"""

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

# ╔═╡ 43101224-7ac5-11eb-104c-0323cf1813c5
md"The zero and one elements are defined as"

# ╔═╡ 792df1aa-7a23-11eb-2991-196336246c43
zero(CountingTropical{Float64})

# ╔═╡ 8388305c-7a23-11eb-1588-79c3c6ce9db9
one(CountingTropical{Float64})

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
		for (edge, J, label) in [((a, b), J12, "J₁₂"), ((b, c), J23, "J₂₃"), ((c, d), J34, "J₃₄"), ((a, d), J14, "J₁₄"), ((d,e), J45, "J₄₅"), ((c, a), J13, "J₁₃")]
			(J > 0 ? eb : eb_) >> edge
		end
	end
	img2 = canvas() do
		for (loc, label) in [(a, "↑"), (b, "↑"), (c, "↓"), (d, "↑"), (e, "↑")]
			sq >> loc
			tb2 >> (loc, label)
		end
		for (edge, J, label) in [((a, b), J12, "J₁₂"), ((b, c), J23, "J₂₃"), ((c, d), J34, "J₃₄"), ((a, d), J14, "J₁₄"), ((d,e), J45, "J₄₅"), ((c, a), J13, "J₁₃")]
			(J > 0 ? eb : eb_) >> edge
		end
	end
	Compose.compose(context(), (context(.05, 0, .5*7/8, 1), img1),
		(context(.55, 0, .5*7/8, 1), img2))
end

# ╔═╡ c2f987aa-7a36-11eb-0d03-4b6d328d8fa4
md"## In our paper"

# ╔═╡ d531c952-7ad9-11eb-1247-dd1913cc4678
html"""<div align="center"><a href="https://github.com/TensorBFS/TropicalTensors.jl">code is available on github <svg class="octicon octicon-mark-github v-align-middle" height="32" viewBox="0 0 16 16" version="1.1" width="32" aria-hidden="true"><path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"></path></svg></a></div>"""

# ╔═╡ 541f4062-7b22-11eb-2eb8-17585a3de9c3
md"#### Square lattices"

# ╔═╡ 080bf23c-7ad8-11eb-37ac-01d2b6439f55
let
	img1 = viz(SquareLattice(32,32); node_style=nodestyle(:circle; r=0.008), text_style=textstyle(:default, fill("transparent")))
	Compose.set_default_graphic_size(14cm*0.4, 7cm*0.8)
	leftright(updown(img1, md"We obtain the exact ground state energy of Ising spin glasses on square lattice up to $32^2$ spins."), updown(HTML("""<img src="https://user-images.githubusercontent.com/6257240/109566189-87bc5980-7ab1-11eb-9d08-99cd573007df.png" width=270px></img>"""), md"""Wall clock time for computing the ground state energy of the (a) Ising spin glass on an open square lattice with
``L^2`` spins. (tensor networks are contracted with [Yao.jl](https://github.com/QuantumBFS/Yao.jl))"""))
end

# ╔═╡ 64f18c2e-7b22-11eb-352f-9d6e228cef49
md"#### Cubic lattices"

# ╔═╡ 9deb5b9a-7adf-11eb-3ba0-0d3716d7d603
let
	Compose.set_default_graphic_size(14cm*0.4, 7cm*0.8)
	θ = 2.31
	ϕ = 2.8
	cam_position = SVector(0.0, 0.0, 0.5)
	rot = RotY(θ)*RotX(ϕ)
	cam_transform = PerspectiveMap() ∘ inv(AffineMap(rot, rot*cam_position))
	Nx = Ny = Nz = 6
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
	leftright(Compose.compose(context(0.5,0.5, 1.0, 1.0), fig), md"""Ground state and degeneracy on cubic lattice up to $6^3$ spins. 
""")
end

# ╔═╡ 72742d66-7b22-11eb-2fac-cb2534558248
md"#### Chimera lattices"

# ╔═╡ 4f3a12e0-7ad5-11eb-2b37-c95342185c3e
let
	Compose.set_default_graphic_size(6cm, 8cm)
	img = Compose.compose(context(0.0,0,8/6,1), viz(ChimeraLattice(8, 8); node_style=nodestyle(:circle; r=0.008), text_style=textstyle(:default, fill("transparent"), fontsize(2pt))))
	leftright(updown(img, md"``\pm J`` Ising spin glass on the chimera graph of D-Wave quantum annealer of $512$ qubits in less than $100$ seconds and investigate the exact value of the residual entropy of $\pm J$ spin glasses on the chimera graph."), updown(html"""<img src="https://user-images.githubusercontent.com/6257240/109566350-bb977f00-7ab1-11eb-953f-127d7919e3e6.png" width=270px/>""", md"Wall clock time for computing the ground
state energy of Ising spin glass on the chimera graph with the
``L \times L`` unit cell (``8L^2`` spins). (tensor networks are contracted with [Yao.jl](https://github.com/QuantumBFS/Yao.jl))"))
end

# ╔═╡ 81e06ac6-7b22-11eb-3042-373a49bbdb49
md"#### Random 3-regular graphs"

# ╔═╡ e59d7a44-7ae7-11eb-3d93-3bc5cc46bc65
let
	function rand_3regular_tn(::Type{T}, n; D=2) where T
		g = LightGraphs.random_regular_graph(n, 3)
		labels = 1:ne(g)
		arrays = [rand(T, fill(2, 3)...) for i=1:n]
		labels = [Int[] for i=1:n]
		for (k, e) in enumerate(LightGraphs.edges(g))
			push!(labels[e.src], k)
			push!(labels[e.dst], k)
		end
		tensors = LabeledTensor.(arrays, labels)
		TensorNetwork(tensors)
	end

	img = SimpleTensorNetworks.viz_tnet(rand_3regular_tn(Float64, 220), node_facecolor="black", linecolor="black", node_fontsize=0)
	#leftright(img, updown(html"""<img src="https://user-images.githubusercontent.com/6257240/109566350-bb977f00-7ab1-11eb-953f-127d7919e3e6.png" width=270px/>""", md"Wall clock time for computing the ground state energy of Ising spin glass on the chimera graph with the ``L \times L`` unit cell (``8L^2`` spins)."))
	Compose.set_default_graphic_size(7cm, 6cm)
	leftright(Compose.compose(context(0.0, 0.0, 6/7, 1.0), img), md"""
The spin glass on the random graphs: our method can compute optimal solutions and count the number of solutions for spin glasses and combinatorial optimization problems on on $3$ regular random graphs up to $220$ spins, on a single GPU. This is inaccessible by existing methods.
""")
end

# ╔═╡ 06bbead0-793f-11eb-0dec-c549b461b9cf
md"""
## Max 2-satisfiability problem
* decision: Polynomial
* counting: #P-complete

A 2-satisfiability problem may be described using a Boolean expression with a special restricted form. It is a conjunction (a Boolean and operation) of clauses, where each clause is a disjunction (a Boolean or operation) of two variables or negated variables. The following is an example of 2-SAT problem of size ``6``

```math
\begin{align}
& (x_{0}\lor x_{2})\land (x_{0}\lor \lnot x_{3})\land (x_{1}\lor \lnot x_{3})\land (x_{1}\lor \lnot x_{4})\land \\
& (x_{2}\lor \lnot x_{4})\land {}(x_{0}\lor \lnot x_{5})\land (x_{1}\lor \lnot x_{5})\land (x_{2}\lor \lnot x_{5})\land \\
&(x_{3}\lor x_{6})\land (x_{4}\lor x_{6})\land (x_{5}\lor x_{6}).
\end{align}
```
(from [wiki](https://en.wikipedia.org/wiki/2-satisfiability))
"""

# ╔═╡ ef2d2446-793f-11eb-223a-c5fe0ed5e367
md"""
We define an energy function as
```math
E = \sum\limits_{k=1}^{|C|} C_k(s_{i_k}, s_{j_k})
```
where ``C_k(s_{i_k}, s_{j_k})`` is ``+1`` if the ``k``th clause on $i$th and $j$th boolean variable is satisfied, otherwise, it is ``-1``.

```math
(T_{e})_{s_i s_j} = C_e(s_i,s_j)
```

```math
T_{v}^n = \delta^{n,q=2}
```
"""

# ╔═╡ 73f517de-7aed-11eb-03d1-db03dfb01a35
md"Since the resulting (counting) tropical number 11 is equal to the number of clauses, all clauses are satisfied, and the degeneracy is 16."

# ╔═╡ 5f2243c4-793d-11eb-1add-392387bb559f
md"""
## Potts model

* decision: NP-Complete
* counting: #P-complete

The Potts model is generalization of Ising model, where a spin takes one of ``q`` possible values, distributed uniformly about the circle, at angles

``\theta _{n}={\frac  {2\pi n}{q}}``,
where ``n = 0, 1, ..., q-1`` and that the interaction Hamiltonian be given by

```math
H=J\sum _{{i,j \in E}}\cos \left(\theta _{{s_{i}}}-\theta _{{s_{j}}}\right)
```
For $q=3$, we have the edge tensor
```math
T_e = J\left(\begin{matrix}1 & -1/2 & -1/2 \\ -1/2 & 1 & -1/2 \\ -1/2 & -1/2 & 1\end{matrix}\right)
```
and the vertex tensor

$$T_v^n=\delta^{n,q=3}$$
"""

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
		for (edge, J, label) in [((a, b), J12, "J₁₂"), ((b, c), J23, "J₂₃"), ((c, d), J34, "J₃₄"), ((a, d), J14, "J₁₄"), ((d,e), J45, "J₄₅"), ((c, a), J13, "J₁₃")]
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

# ╔═╡ 5f4e0fec-7afd-11eb-37c7-11b84027136a
md"There are 4 possible configurations"

# ╔═╡ f13469bc-7afb-11eb-3dab-2b6cdf290f6f
let
	sq = nodestyle(:circle, fill("black"); r=0.04)
	sq_ = nodestyle(:circle, fill("white"), stroke("black"); r=0.04)
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
	nodes = [a,b,c,d,e]
	function drawconfig(config)
		canvas() do
			for (loc, c) in zip(nodes, config)
				(c==1 ? sq : sq_) >> loc
			end
			for (edge, label) in [((a, b), "+1"), ((b, c), "-1"), ((c, d), "-1"), ((a, d), "+1"), ((d,e), "+1"), ((c, a), "+1")]
				eb >> edge
			end
		end
	end
	Compose.set_default_graphic_size(16cm, 4cm)
	img1 = drawconfig([1, 0, 0, 0, 1])
	img2 = drawconfig([0, 1, 0, 0, 1])
	img3 = drawconfig([0, 0, 1, 0, 1])
	img4 = drawconfig([0, 1, 0, 1, 0])
	Compose.compose(context(),
	(context(0, 0, 0.25, 1), img1),
	(context(0.25, 0, 0.25, 1), img2),
	(context(0.5, 0, 0.25, 1), img3),
	(context(0.75, 0, 0.25, 1), img4),
	)
end

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

# ╔═╡ 5caa1b7c-7b21-11eb-394d-379351fe5170
md"Then one can also get the correct result with einsum."

# ╔═╡ 1b08b3ac-7b1e-11eb-2249-ddd787c549d4
let
	T = CountingTropical{Float64}
	ein"ab,bc,cd,ad,de,ac->"(
		T.([0 1; 1 -Inf]),  # (1, 2)
		T.([0 1; 0 -Inf]),	# (2, 3)
		T.([0 1; 0 -Inf]),	# (3, 4)
		T.([0 0; 0 -Inf]),	# (1, 4)
		T.([0 1; 0 -Inf]),	# (4, 5)
		T.([0 0; 0 -Inf]),	# (1, 3)
	)[]
end

# ╔═╡ 7bdf517e-79ff-11eb-38a3-49c02d94d943
md"## The Song Shan Lake Spring School (SSSS) Challendge"

# ╔═╡ a387b018-79ff-11eb-0383-1fcf82853afc
md"""
[You may find more than 10 solutions in our github repo](https://github.com/QuantumBFS/SSSS/blob/master/Challenge.md)
"""

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
end;

# ╔═╡ 1dbb9e90-78b0-11eb-2014-6dc6cfb35387
md"``R_x`` = $(@bind θ2 Slider(0:0.01:2π; default=0.5, show_value=true))"

# ╔═╡ 1dbc9afc-78b0-11eb-0940-2dcadf5408bb
md"``R_y`` = $(@bind ϕ2 Slider(0:0.01:2π; default=2.8, show_value=true))"

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

# ╔═╡ 88e14ef2-7af1-11eb-23d6-b34b1eff8f87
md"""
In the Buckyball structure as shown in the figure, we attach an ising spin ``s_i=\pm 1`` on each vertex. The neighboring spins interact with an anti-ferromagnetic coupling of unit strength.

1. Get ``\ln Z/N``, where ``N`` is the number of vertices, and
```math
Z = \sum_{\{s\}}e^{-\sum_{i,j \in E}s_i s_j}
```

2. Count the ground state degeneracy.
"""

# ╔═╡ b6560404-7b2d-11eb-21d7-a1e55609ebf7
# the positions of fullerene atoms
c60_xy = fullerene();

# ╔═╡ 6f649efc-7b2d-11eb-1e80-53d84ef98c13
# find edges: vertex pairs with distance smaller than 5.
c60_edges = [(i=>j) for (i,(i2,j2,k2)) in enumerate(c60_xy), (j,(i1,j1,k1)) in enumerate(c60_xy) if i<j && (i2-i1)^2+(j2-j1)^2+(k2-k1)^2 < 5.0];

# ╔═╡ 20125640-79fd-11eb-1715-1d071cc6cf6c
md"construct tensor network by assigning labels `(edge index, boolean)` to tensors, where the boolean identifies whether this label correspond to the source node or the destination side of the edge. The resulting tensor network contains 90 edge tensors and 60 vertex tensors."

# ╔═╡ 698a6dd0-7a0e-11eb-2766-1f0baa1317d2
md"find a proper contraction order by greedy search"

# ╔═╡ 1c4b19d2-7b30-11eb-007b-ab03052b22d2
md"The greedy contraction order can be visualized by dragging the slider (if you run it on your local host)"

# ╔═╡ 4c137484-7b30-11eb-2fb1-190d8beebbc3
md"Finding the optimal contraction order is know as NP-hard, however, there are [some heuristic algorithms](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.060503) can find good contraction path in limited time.
"

# ╔═╡ e302bd1c-7ab5-11eb-03f6-69dcbb817354
md"# Future Directions

* Finding the optimal contracting order of `einsum`,
    - einsum has less redundancy, with less effort in mapping, and potentially faster in contraction,
    - tensor network can utilize BLAS easily,
    - the contraction order for tensor networks are more studies.
* TropicalBLAS,
    - TropicalGEMM on CPU for Tropical numbers (DONE),
    - TropicalGEMM on GPU (WIP),
    - TropicalGEMM on CPU for CountingTropical numbers,
"

# ╔═╡ d53a5e5a-7b33-11eb-0a65-4d4844e9cf0a
md"# Utilities"

# ╔═╡ e51dc3e8-7b33-11eb-12ee-57413f6b8bca
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

# ╔═╡ f914a760-7ad2-11eb-17e6-c39cf676196e
# a `δ` tensor with rank 2 and size 3 × 3
δtensor(Tropical{Float64}, 3, 2)

# ╔═╡ 37472f2a-7a2a-11eb-1be3-13513d61fcb2
# we use the `@ein_str` macro from OMEinsum to perform the contraction
# Here we use a, b, c, d, e to represent the original vertex labels 1, 2, 3, 4, 5
ein"ab,bc,cd,ad,de,ac->"([ising_bondtensor(Tropical{Float64}, J) for J in [J12, J23, J34, J14, J45, J13]]...)[]

# ╔═╡ 15cf0c36-7a21-11eb-3e14-63950bcce943
tnet = let
	T = Tropical{Float64}
	edges = [(1, 2), (2, 3), (3, 4), (1, 4), (4, 5), (1, 3)]
	build_tensornetwork(
		vertices = 1:5,
		vertex_arrays = [ising_vertextensor(T, count(e->i ∈ e, edges), 0.0) for i=1:5],
		edges = edges,
		edge_arrays = [ising_bondtensor(T, J) for J in [J12, J23, J34, J14, J45, J13]]
	)
end;

# ╔═╡ a9544e66-7a27-11eb-2b27-1d2124988fb2
contraction_result = let
	tc, sc, trees = trees_greedy(tnet; strategy="min_reduce")
	SimpleTensorNetworks.contract(tnet, trees[1]).array[]
end

# ╔═╡ 3bb2e0c2-7a28-11eb-1ea5-ab03d16bf0b3
md"The mininum energy is $(-Int(contraction_result.n))."

# ╔═╡ c1a7bd4a-7a36-11eb-176a-f399eb6b5f49
# `CountingTropical{Float64}` has two fields `n` (number) and `c` (counting), both are of type `Float64`.
ein"ab,bc,cd,ad,de,ac->"([ising_bondtensor(CountingTropical{Float64}, J) for J in [J12, J23, J34, J14, J45, J13]]...)[]

# ╔═╡ b8c8999a-7aec-11eb-3ccd-69b48fcb93c2
let
	T = CountingTropical{Float64}
	ein"ac,ad,bd,be,ce,af,bf,cf,dg,eg,fg->"(
		twosat_bondtensor(T, true, true),
		twosat_bondtensor(T, true, false),
		twosat_bondtensor(T, true, false),
		twosat_bondtensor(T, true, false),
		twosat_bondtensor(T, true, false),
		twosat_bondtensor(T, true, false),
		twosat_bondtensor(T, true, false),
		twosat_bondtensor(T, true, false),
		twosat_bondtensor(T, true, true),
		twosat_bondtensor(T, true, true),
		twosat_bondtensor(T, true, true),
	)[]
end

# ╔═╡ 0405f4d8-7afb-11eb-2163-597b2edcf17e
tensor_network_mis = let
	T = CountingTropical{Float64}
	edges = [(1, 2), (2, 3), (3, 4), (1, 4), (4, 5), (1, 3)]
	build_tensornetwork(
		vertices = 1:5,
		vertex_arrays = [mis_vertextensor(T, count(e->i ∈ e, edges)) for i=1:5],
		edges = edges,
		edge_arrays = [mis_bondtensor(T) for i=1:6]
	)
end;

# ╔═╡ 3a34aaec-7afb-11eb-1fc4-2fbc027753cc
contraction_result_mis = let
	tc, sc, trees = trees_greedy(tensor_network_mis; strategy="min_reduce")
	SimpleTensorNetworks.contract(tensor_network_mis, trees[1]).array[]
end

# ╔═╡ c26b5bb6-7984-11eb-18fe-2b6a524f5c85
c60_tnet = let
	T = CountingTropical{Float64}
	build_tensornetwork(
		vertices=1:60,
		vertex_arrays = [ising_vertextensor(T, 3, 0.0) for j=1:length(c60_xy)],
		edges = c60_edges,
		edge_arrays = [ising_bondtensor(T, -1.0) for i = 1:length(c60_edges)]
	)
end;

# ╔═╡ ae92d828-7984-11eb-31c8-8b3f9a071c24
tcs, scs, c60_trees = (Random.seed!(2); trees_greedy(c60_tnet; strategy="min_reduce"));

# ╔═╡ 2b899624-798c-11eb-20c4-fd5523f7abff
md"time complexity = $(round(log2sumexp2(tcs); sigdigits=4)), space complexity = $(round(maximum(scs); sigdigits=4))"

# ╔═╡ d2161642-798a-11eb-2dec-cfe6cda6af5c
SimpleTensorNetworks.contract(c60_tnet, c60_trees[]).array[]

# ╔═╡ 58e38656-7b2e-11eb-3c70-25a919f9926a
md"contraction step = $(@bind nstep_c60 Slider(0:length(c60_tnet); show_value=true, default=60))"

# ╔═╡ 12740186-7b2f-11eb-35e4-01e6f9ffbb4d
c60_contraction_masks = let
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
	contraction_mask(c60_tnet, c60_trees[])
end;

# ╔═╡ c1c74e70-7b2c-11eb-2f26-21f54ad00fb2
let
	mask = c60_contraction_masks[nstep_c60+1]
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
	Compose.compose(context(0.5,0.35, 1.0, 1.0), fig)
end

# ╔═╡ Cell order:
# ╟─c456b902-7959-11eb-03ba-dd14a2cd5758
# ╟─dfa8834c-e8c6-49b4-8bde-0816b573cbee
# ╟─121b4926-7aba-11eb-30e1-7b8edd4f0166
# ╠═5bb40ad6-7b33-11eb-0b31-63d5e47fa0e7
# ╟─3205a536-7a17-11eb-3473-b71305c96ca4
# ╟─3208fd8a-7a17-11eb-35ce-4d6b141c1aff
# ╟─32116a92-7a17-11eb-228f-0713510d0348
# ╟─3217d1ca-7a17-11eb-14eb-a77ccfd983a0
# ╟─3221a326-7a17-11eb-0fe6-f75798a411b9
# ╟─32277c3a-7a17-11eb-3763-af68dbb81465
# ╟─322d2958-7a17-11eb-2deb-613b7680a5bb
# ╟─3237e33e-7a17-11eb-2869-b92d0801bc6e
# ╟─ec841be8-7a16-11eb-3337-376e26b7da25
# ╟─be76e52a-7852-11eb-179b-afbc6efcab55
# ╟─d0b54b76-7852-11eb-2398-0911380fa090
# ╠═af13e090-7852-11eb-21ae-8b94f25f1a4f
# ╠═d770f232-7864-11eb-0e9a-81528e359d39
# ╟─5d16a472-785b-11eb-1b94-dd6d8f860c24
# ╟─7da42d08-7ad2-11eb-095a-87aedba91b35
# ╟─3372871c-785b-11eb-3092-4bbc419cb788
# ╠═2173a6cc-785b-11eb-1ab6-7fb875224dd9
# ╟─518b7d4e-785b-11eb-3b7c-1389065b9cbd
# ╠═2868b292-785b-11eb-015e-6b5613bd9e39
# ╟─98ae0960-797d-11eb-3646-c5b7e05d3f7c
# ╠═f914a760-7ad2-11eb-17e6-c39cf676196e
# ╟─86921d00-7a17-11eb-2695-add5f9eeda5b
# ╟─915e8096-7a17-11eb-177d-a39ffed7ca91
# ╟─958c489e-7a17-11eb-2be4-cfdb44da7d2f
# ╟─3d7ca08c-7b01-11eb-1d78-af35dc7e577c
# ╟─368c9728-7b01-11eb-3c30-c9d8b4d00ace
# ╟─25d64fd4-7b04-11eb-105a-71f98c236ec6
# ╟─0e891c90-7b0d-11eb-2616-83f1e9a21ae2
# ╠═71ae3014-7b01-11eb-1707-2f6d249d8cf6
# ╟─14b93998-7b25-11eb-22cc-f1a34c46570f
# ╟─48749c00-7b25-11eb-19f6-41b57a9808ab
# ╠═8bd6eee8-7b03-11eb-2a79-694e7bc29bd6
# ╟─5da00572-7b04-11eb-01f1-cd94579b478e
# ╟─211911da-7a18-11eb-12d4-65b0dec4b8dc
# ╟─281a5cf0-7a18-11eb-3385-c3e64f41e4da
# ╟─2c47b692-7a18-11eb-2878-2976435507f5
# ╟─442bcb3c-7940-11eb-18e5-d3158b74b1dc
# ╟─f7208b6e-793c-11eb-0dfa-0d63752ba53e
# ╟─d3b4f162-7ad4-11eb-271c-677cb307c447
# ╟─22118a36-7a36-11eb-18c3-dd2adac6118b
# ╠═c85217b6-7a23-11eb-04c5-fb4dc9c02ef1
# ╟─b3b986aa-7a1d-11eb-17d4-e5675015b221
# ╟─b5cd769e-7a1e-11eb-1d82-e1c265dfdd52
# ╟─00ed185e-7a2d-11eb-1b27-cb834e75e916
# ╟─9e6fbf32-7a2f-11eb-17cb-9167d6a34281
# ╟─6e1c507e-7a1a-11eb-05bc-dbccc3aebdf9
# ╟─b52ead96-7a2a-11eb-334f-e5e5ff5867e3
# ╟─05109d30-7a29-11eb-320a-fb0b0d8e2632
# ╟─d0ecd3f2-7a2d-11eb-126d-7dab740d8e1f
# ╟─9c860e2a-7a2e-11eb-231f-63e9aca1daa0
# ╟─5efee244-7a2d-11eb-3782-b9d55086d623
# ╠═1749c0f2-7a2a-11eb-1932-07a7f920b0da
# ╠═37472f2a-7a2a-11eb-1be3-13513d61fcb2
# ╟─023ebf7c-7b36-11eb-1c9f-430773395534
# ╟─96290770-7a20-11eb-0ac8-33a6492c7b12
# ╟─c1f90d6c-7a1d-11eb-2843-f971b5f6f3b0
# ╟─64e08a56-7a36-11eb-29fd-03662b4d6612
# ╟─f54119ca-7a1e-11eb-1bec-bf855e34658d
# ╟─37102544-7abf-11eb-3ac4-6702dfc55425
# ╠═15cf0c36-7a21-11eb-3e14-63950bcce943
# ╠═a9544e66-7a27-11eb-2b27-1d2124988fb2
# ╟─3bb2e0c2-7a28-11eb-1ea5-ab03d16bf0b3
# ╟─695e405c-786d-11eb-0a6e-bb776d9626ad
# ╟─43101224-7ac5-11eb-104c-0323cf1813c5
# ╠═792df1aa-7a23-11eb-2991-196336246c43
# ╠═8388305c-7a23-11eb-1588-79c3c6ce9db9
# ╠═c1a7bd4a-7a36-11eb-176a-f399eb6b5f49
# ╟─4190393a-7ac4-11eb-3ac6-eb8e3574fdc9
# ╟─56fdb22c-7ac4-11eb-2831-a777d9ca89f3
# ╟─c2f987aa-7a36-11eb-0d03-4b6d328d8fa4
# ╟─d531c952-7ad9-11eb-1247-dd1913cc4678
# ╟─541f4062-7b22-11eb-2eb8-17585a3de9c3
# ╟─080bf23c-7ad8-11eb-37ac-01d2b6439f55
# ╟─64f18c2e-7b22-11eb-352f-9d6e228cef49
# ╟─9deb5b9a-7adf-11eb-3ba0-0d3716d7d603
# ╟─72742d66-7b22-11eb-2fac-cb2534558248
# ╟─4f3a12e0-7ad5-11eb-2b37-c95342185c3e
# ╟─81e06ac6-7b22-11eb-3042-373a49bbdb49
# ╟─e59d7a44-7ae7-11eb-3d93-3bc5cc46bc65
# ╟─06bbead0-793f-11eb-0dec-c549b461b9cf
# ╟─ef2d2446-793f-11eb-223a-c5fe0ed5e367
# ╠═b8c8999a-7aec-11eb-3ccd-69b48fcb93c2
# ╟─73f517de-7aed-11eb-03d1-db03dfb01a35
# ╟─5f2243c4-793d-11eb-1add-392387bb559f
# ╟─e739e74c-7af0-11eb-104f-5f94da1bf0be
# ╟─80d6c2b6-7aef-11eb-1bf5-5d4f266dfa73
# ╟─344042b4-793d-11eb-3d6f-43eb2a4db9f4
# ╟─80d764a8-7afd-11eb-3fb8-79169ca56c7e
# ╟─04d11828-7afa-11eb-3e73-1bbecf566f74
# ╟─5d95a598-7afa-11eb-10eb-db79fa44dd2a
# ╟─d29470d4-7afa-11eb-0afc-a34e39d49aa5
# ╟─75c37046-7b1b-11eb-00f5-7fc49f73f4d9
# ╠═0405f4d8-7afb-11eb-2163-597b2edcf17e
# ╠═3a34aaec-7afb-11eb-1fc4-2fbc027753cc
# ╟─5f4e0fec-7afd-11eb-37c7-11b84027136a
# ╟─f13469bc-7afb-11eb-3dab-2b6cdf290f6f
# ╟─891c39b8-7b1b-11eb-2d70-67fd24021027
# ╟─5caa1b7c-7b21-11eb-394d-379351fe5170
# ╠═1b08b3ac-7b1e-11eb-2249-ddd787c549d4
# ╟─7bdf517e-79ff-11eb-38a3-49c02d94d943
# ╟─a387b018-79ff-11eb-0383-1fcf82853afc
# ╟─5a5d4de6-7895-11eb-15c6-bda7a4342002
# ╟─1dbb9e90-78b0-11eb-2014-6dc6cfb35387
# ╟─1dbc9afc-78b0-11eb-0940-2dcadf5408bb
# ╟─9b1dc21a-7896-11eb-21f6-bfe9b4dc9ccf
# ╟─88e14ef2-7af1-11eb-23d6-b34b1eff8f87
# ╠═b6560404-7b2d-11eb-21d7-a1e55609ebf7
# ╠═6f649efc-7b2d-11eb-1e80-53d84ef98c13
# ╟─20125640-79fd-11eb-1715-1d071cc6cf6c
# ╠═c26b5bb6-7984-11eb-18fe-2b6a524f5c85
# ╟─698a6dd0-7a0e-11eb-2766-1f0baa1317d2
# ╠═ae92d828-7984-11eb-31c8-8b3f9a071c24
# ╟─2b899624-798c-11eb-20c4-fd5523f7abff
# ╠═d2161642-798a-11eb-2dec-cfe6cda6af5c
# ╟─1c4b19d2-7b30-11eb-007b-ab03052b22d2
# ╟─58e38656-7b2e-11eb-3c70-25a919f9926a
# ╟─12740186-7b2f-11eb-35e4-01e6f9ffbb4d
# ╟─c1c74e70-7b2c-11eb-2f26-21f54ad00fb2
# ╟─4c137484-7b30-11eb-2fb1-190d8beebbc3
# ╟─e302bd1c-7ab5-11eb-03f6-69dcbb817354
# ╟─d53a5e5a-7b33-11eb-0a65-4d4844e9cf0a
# ╠═e51dc3e8-7b33-11eb-12ee-57413f6b8bca
