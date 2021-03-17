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
	
	function viz_lattice(lt::Viznet.AbstractSites; ncolor, r)
		line_style=bondstyle(:default, stroke("#333333"), linewidth(4*unit(lt)*mm))
		colors = if ncolor==1
			["#333333"]
		elseif ncolor==2
			["#CC4040", "#4040CC"]
		elseif ncolor==3
			["#CC4040", "#40CC40", "#4040CC"]
		else
			error("!!")
		end
		node_styles = [nodestyle(:circle, fill(color), r=r) for color in colors]
		#text_style=textstyle(:default, fontsize(unit(lt)*100pt)), labels=vertices(lt)
		canvas() do
			for (i,node) in enumerate(Viznet.vertices(lt))
				rand(node_styles) >> lt[node]
				#text_style >> (lt[node], "$(labels[i])")
			end
			for bond in Viznet.bonds(lt)
				line_style >> lt[bond[1]; bond[2]]
			end
		end
	end
	
	
	function _spring_layout(adj_matrix::AbstractMatrix,
						   locs_x=2*rand(size(adj_matrix, 1)).-1.0,
						   locs_y=2*rand(size(adj_matrix, 1)).-1.0;
						   C=2.0,
						   MAXITER=100,
						   INITTEMP=2.0)

		nvg = size(adj_matrix, 1)

		# The optimal distance bewteen vertices
		k = C * sqrt(4.0 / nvg)
		k² = k * k

		# Store forces and apply at end of iteration all at once
		force_x = zeros(nvg)
		force_y = zeros(nvg)

		# Iterate MAXITER times
		@inbounds for iter = 1:MAXITER
			# Calculate forces
			for i = 1:nvg
				force_vec_x = 0.0
				force_vec_y = 0.0
				for j = 1:nvg
					i == j && continue
					d_x = locs_x[j] - locs_x[i]
					d_y = locs_y[j] - locs_y[i]
					dist²  = (d_x * d_x) + (d_y * d_y)
					dist = sqrt(dist²)

					if !( iszero(adj_matrix[i,j]) && iszero(adj_matrix[j,i]) )
						# Attractive + repulsive force
						# F_d = dist² / k - k² / dist # original FR algorithm
						F_d = dist / k - k² / dist²
					else
						# Just repulsive
						# F_d = -k² / dist  # original FR algorithm
						F_d = -k² / dist²
					end
					force_vec_x += F_d*d_x
					force_vec_y += F_d*d_y
				end
				force_x[i] = force_vec_x
				force_y[i] = force_vec_y
			end
			# Cool down
			temp = INITTEMP / iter
			# Now apply them, but limit to temperature
			for i = 1:nvg
				fx = force_x[i]
				fy = force_y[i]
				force_mag  = sqrt((fx * fx) + (fy * fy))
				scale      = min(force_mag, temp) / force_mag
				locs_x[i] += force_x[i] * scale
				locs_y[i] += force_y[i] * scale
			end
		end

		# Scale to unit square
		min_x, max_x = minimum(locs_x), maximum(locs_x)
		min_y, max_y = minimum(locs_y), maximum(locs_y)
		function scaler(z, a, b)
			2.0*((z - a)/(b - a)) - 1.0
		end
		map!(z -> scaler(z, min_x, max_x), locs_x, locs_x)
		map!(z -> scaler(z, min_y, max_y), locs_y, locs_y)

		return [((x+1)/2, (y+1)/2) for (x, y) in zip(locs_x, locs_y)]
	end

	function viz_graph(graph; r=0.25/sqrt(nv(graph)+1), show_edgeindex=false,
        node_fontsize=100pt/sqrt(nv(graph)+1),
        edge_fontsize=200pt/sqrt(nv(graph)+1),
        labels=1:nv(graph),
        locs=_spring_layout(adjacency_matrix(graph)),
        linecolor="#333333",
        node_edgecolor="transparent",
        node_facecolors=["#CC4040", "#4040CC"]
    )
		nt = nv(graph)
		nbs = [nodestyle(:default, fill(c), stroke(node_edgecolor), linewidth(2mm/sqrt(nv(graph)+1)); r=r) for c in node_facecolors]
		eb = bondstyle(:default, linewidth(3mm/sqrt(nv(graph)+1)), stroke(linecolor))
		Compose.compose(Compose.context(r, r, 1-2r, 1-2r), canvas() do
			for (loc, label) in zip(locs, labels)
				rand(nbs) >> loc
			end
			for edge in edges(graph)
				eb >> (locs[edge.src], locs[edge.dst])
			end
		end)
	end


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
using TropicalNumbers,  		# tropical number type
		LightGraphs,			# graph operations
		Random,
		OMEinsum,				# Einstein's summation notation
    	SimpleTensorNetworks  	# tensor network contraction

# ╔═╡ 3674a622-823b-11eb-3991-d5771010237b
using Pkg; Pkg.status()

# ╔═╡ 94b870d2-8235-11eb-33e7-35bf5132efd6
using Profile

# ╔═╡ dfa8834c-e8c6-49b4-8bde-0816b573cbee
html"""
<style>
body {
counter-reset: section subsection example}

h2::before {
counter-reset: subsection;
  counter-increment: section;
  content: counter(section) ". ";
}
</style>
<div align="center">
<a class="Header-link " href="https://github.com/TensorBFS/TropicalTensors.jl" data-hotkey="g d" aria-label="Homepage " data-ga-click="Header, go to dashboard, icon:logo">
  <svg class="octicon octicon-mark-github v-align-middle" height="32" viewBox="0 0 16 16" version="1.1" width="32" aria-hidden="true"><path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"></path></svg>
</a>
<br>
<a href="https://raw.githubusercontent.com/GiggleLiu/notebooks/master/notebooks/tropical/tropicaltensornetwork.jl" target="_blank"> download this notebook </a></div>
"""

# ╔═╡ 121b4926-7aba-11eb-30e1-7b8edd4f0166
html"""<h1>Tropical tensor networks</h1>
<p><big>for solving combinatoric optimization problems efficiently</big></p>

<a href='https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.126.090506'>Phys. Rev. Lett. 126, 090506 (2021)</a></p>
<p>Jin-Guo Liu, Lei Wang, and Pan Zhang</p>
"""

# ╔═╡ 3205a536-7a17-11eb-3473-b71305c96ca4
md"## What is a tropical tensor network?"

# ╔═╡ 592825a5-890a-4de8-984f-3d4ca48daca4
md"A Tropical tensor network is a tensor network with Tropical numbers inside."

# ╔═╡ 5066553a-9b22-4f3b-b8ba-13a88291d9b9
md"**What is a tensor network?**"

# ╔═╡ 3208fd8a-7a17-11eb-35ce-4d6b141c1aff
md"Tensor network is a generalization of matrix multiplication to multiple tensors. The graphical representation of matrix multiplication is
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
* a matrix is a tensor of rank-2, it is represented as a vertex with two labeled edges (or legs),
* the **contraction** (a generalization of matrix multiplication to tensors) is represented by connecting legs with same labels, where connecting a pair of legs means **sum-product** over a specific degree of freedom corresponding to the leg label,
"

# ╔═╡ 3221a326-7a17-11eb-0fe6-f75798a411b9
md"We replace matrices with tensors, matrix multiplication with tensor contraction and allow more than two nodes in the graph, and we will get a tensor network. The following is an example of graphical representation of a tensor network
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
*To know more about tensor networks:*

* A Practical Introduction to Tensor Networks: Matrix Product States and Projected Entangled Pair States, [arXiv: 1306.2164](https://arxiv.org/abs/1306.2164)"

# ╔═╡ ec841be8-7a16-11eb-3337-376e26b7da25
md"**What are tropical numbers?**"

# ╔═╡ d0b54b76-7852-11eb-2398-0911380fa090
md"""
Tropical numbers are numbers with tropical algebra. Tropical algebra is defined by replacing the usual sum and product operators for ordinary real numbers with the max and sum operators respectively 

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
**Tropical algebra is a semi-ring!**

The $\oplus$ and $\odot$ operators still have commutative, associative, and distributive properties. However, since there is no additive inverse, the $\oplus$ and $\odot$ and operations define a semiring over ``\mathbb R \cup \{-\infty\}``. 
"""

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
The tropical $\delta$ tensor (the generalization of identity matrix to tensors) of rank $n$ and dimension $q$ is defined as
```math
δ_{s_i s_j\ldots s_n}^{n,q} = \begin{cases}
 0, & s_i = s_j =\ldots s_n\\
 -\infty, &otherwise
\end{cases}
```
where $s_i,s_j,\ldots s_n \in \{1,2,\ldots q\}$.
"""

# ╔═╡ ca9a13b8-0272-4557-97d4-9168923d363f
# The first parameter is for specifying the element type. It is a common design pattern in Julia, please check https://docs.julialang.org/en/v1/manual/methods/#Extracting-the-type-parameter-from-a-super-type
function δtensor(::Type{T}, q::Int, n::Int) where {T}
	res = zeros(T, fill(q, n)...)
	for i=1:q
		res[fill(i, n)...] = one(T)
	end
	res
end

# ╔═╡ f914a760-7ad2-11eb-17e6-c39cf676196e
# a `δ` tensor with rank 2 and size 3 × 3
δtensor(Tropical{Float64}, 3, 2)

# ╔═╡ 86921d00-7a17-11eb-2695-add5f9eeda5b
md"**Tropical matrix multiplication**

By putting tropical numbers inside matrices, we get tropical matrix multiplication.
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
md"Tropical matrix multiplication is directly related to finding the shortest path. Consider the following graph"

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
*To know more about applications of tropical matrices*

* [Tropical Arithmetic and Shortest Paths](https://personalpages.manchester.ac.uk/staff/mark.muldoon/Teaching/DiscreteMaths/LectureNotes/TropicalShortestPaths.pdf)

* [Methods and Applications of (max,+) Linear Algebra](https://link.springer.com/chapter/10.1007/BFb0023465)"

# ╔═╡ 211911da-7a18-11eb-12d4-65b0dec4b8dc
md"
**Tropical tensor networks**

By putting Tropical numbers inside tensor networks, we get
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
<h2>Mapping hard problems to Tropical Tensor networks</h2>
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
md"""### Ising Spin glass
* decision: NP-Complete
* counting: #P-complete

Ising spin glass is an energy model defined on a graph. It is about finding the maximum of the following energy function
```math
-E = \sum\limits_{i,j\in E} J_{ij} s_i s_j + \sum\limits_{i\in V} h_{i} s_i
```
where $s_i,s_j,\ldots s_n \in \{-1,1\}$, this problem is hard because there are $2^n$ possible configurations of $\{s_1,s_2\ldots s_n\}$. We define ``-E`` rather than ``E`` directly because people are more interested to know the ground state (state with the lowest energy) rather than the highest energy.
The decision version of Ising spin glass is NP-complete, which mean it is unlikely to find a solution to the problem *"what is the lowest possible energy?"* in polynomial time on a classical computer. Its counting version asks *"what is the number of possible configurations that gives the lowest energy"*. It is #P-complete, which is at least as difficult as its dicision version.




"""

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
md"The goal is to find an optimal assignment of $s_1, s_2,\ldots s_n$ that minimizes the energy. Mapping the problem to Einsum is straight-forward. We can rewrite the energy definition in tropical contraction format"

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

# ╔═╡ 624f57db-7f07-4281-a547-d229b9a8413a
function ising_bondtensor(::Type{T}, J) where T
	e = T(J)
	e_ = T(-J)
	[e e_; e_ e]
end

# ╔═╡ 05109d30-7a29-11eb-320a-fb0b0d8e2632
md"where a spin $s_i = 1$ (or $s_i=-1$) is equivalent to $1$ (or $2$) when used in matrix indexing. One can easily check it is the same as the previous negative energy function. The graphical representation is"

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

# ╔═╡ 37472f2a-7a2a-11eb-1be3-13513d61fcb2
# we use the `@ein_str` macro from OMEinsum to perform the contraction
# Here we use a, b, c, d, e to represent the original vertex labels 1, 2, 3, 4, 5
ein"ab,bc,cd,ad,de,ac->"([ising_bondtensor(Tropical{Float64}, J) for J in [J12, J23, J34, J14, J45, J13]]...)[]

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

# ╔═╡ b975680f-0b78-4178-861f-5da6d10327e4
function ising_vertextensor(::Type{T}, n::Int, h) where T
	res = zeros(T, fill(2, n)...)
	res[1] = T(h)
	res[end] = T(-h)
	return res
end

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
html"""<p>Small circles and big circles are for vertex tensors and edge tensors respectively.</p>"""

# ╔═╡ 25553aa4-eef8-40d4-bcf6-0a27f50ab1da
md"""Then we construct tensor network by assigning labels `(edge index, boolean)` ($(HTML("<span style='color:green'>green</span>")) texts in the graph) to tensors, where the boolean identifies whether this label correspond to the source node or the destination side of the edge. """

# ╔═╡ 75764e5b-d2f2-4178-a00f-8f362b59a0b8
function build_tensornetwork(; vertices, vertex_arrays, edges, edge_arrays)
	TensorNetwork([
	# vertex tensors
	[LabeledTensor(vertex_arrays[i], [(j, v==e[1]) for (j, e) in enumerate(edges) if v ∈ e]) for (i, v) in enumerate(vertices)]...,
	# bond tensors
	[LabeledTensor(edge_arrays[j], [(j, true), (j, false)]) for j=1:length(edges)]...
])
end

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

# ╔═╡ a4bf7de3-7794-4609-88e7-9092f496a7bb
md"The tensor network is then contracted with the greedy algorithm."

# ╔═╡ a9544e66-7a27-11eb-2b27-1d2124988fb2
contraction_result = let
	# trees_greedy returns a tuple of (time complexity, space complexity, contraction tree)
	tc, sc, trees = trees_greedy(tnet; strategy="min_reduce")
	SimpleTensorNetworks.contract(tnet, trees[1]).array[]
end

# ╔═╡ 3bb2e0c2-7a28-11eb-1ea5-ab03d16bf0b3
md"The mininum energy is $(-Int(contraction_result.n))."

# ╔═╡ 695e405c-786d-11eb-0a6e-bb776d9626ad
md"
##### Using Tropical numbers for counting
Moreover, one can also employ the present approach to count the number of ground states at the same computational complexity of computing the ground state energy. To implement this, we further generalize the tensor element to be a tuple $(n, c)$ composed by a tropical number $x$ and an ordinary number $n$. The tropical number records the negative energy, while the ordinary number counts the number of minimal energy configurations. For tensor network contraction, we need the multiplication and addition of the tuple: $(n_1, c_1) \odot (n_2,c_2) = (n_1 + n_2, c_1\cdot c_2)$ and 
```math
\begin{equation}
    (n_1, c_1)\oplus (n_2, c_2) = \begin{cases}
 (n_1\oplus n_2, \, c_1 + c_2 ) & \text{if $n_1 = n_2$} \\
 (n_1\oplus n_2,\, c_1 ) & \text{if $n_1>n_2$} \\
 (n_1\oplus n_2,\, c_2 )& \text{if $n_1 < n_2$}
 \end{cases}.
\end{equation}
```
"

# ╔═╡ 1bb36c52-a171-4993-ac86-2250e1e87a01
md"It corresponds to the following four processes of concatenating and comparing configrations on graphs (or tensor networks)."

# ╔═╡ ade34905-5a61-4f71-b347-e02fab120b5d
let
	Compose.set_default_graphic_size(15cm, 10cm)
	a = 0.1
	b = 0.05
	nodes = [(-a, -b), (a, -b), (a, b), (-a, b)]
	nb = Compose.compose(context(), polygon(nodes), stroke("black"), fill("white"))
	tb = textstyle(:default)
	tt = Compose.compose(context(), text(0.0, 0.0, ""))
	x_title = 0.1
	canvas() do
		for (y_1, op, title) in zip(
				[0.2, 0.4, 0.6, 0.8],
				["⊙", "⊕", "⊕", "⊕"],
				[
				("concatenate best configurations of two subgraphs", "value = n₁\ndegeneracy = c₁", "value = n₂\ndegeneracy = c₂", "value = n₁ + n₂\ndegeneracy = c₁ * c₂"),
				("compare two configurations: case n₁ == n₂", "value = n\ndegeneracy = c₁", "value = n\ndegeneracy = c₂", "value = n\ndegeneracy = c₁ + c₂"),
				("compare two configurations: case n₁ > n₂", "value = n₁\ndegeneracy = c₁", "value = n₂\ndegeneracy = c₂", "value = n₁ \ndegeneracy = c₁"),
				("compare two configurations: case n₁ < n₂", "value = n₁\ndegeneracy = c₁", "value = n₂\ndegeneracy = c₂", "value = n₂\ndegeneracy = c₂")
				]
			)
			x = (0.2, y_1)
			y = (0.5, y_1)
			z = (0.8, y_1)
			nb >> x
			nb >> y
			nb >> z
			tt >> ((x_title, y_1-0.08), title[1])
			tb >> (x, title[2])
			tb >> (y, title[3])
			tb >> ((x .+ y) ./ 2, op)
			tb >> ((z .+ y) ./ 2, "=")
			tb >> (z, title[4])
		end
	end
end

# ╔═╡ 43101224-7ac5-11eb-104c-0323cf1813c5
md"The zero and one elements are defined as"

# ╔═╡ 792df1aa-7a23-11eb-2991-196336246c43
zero(CountingTropical{Float64})

# ╔═╡ 8388305c-7a23-11eb-1588-79c3c6ce9db9
one(CountingTropical{Float64})

# ╔═╡ 6faf5bf0-a8f5-4db0-839a-9389d541a690
md"By replacing the tropical numbers with the counting tropical numbers in our privious program, we get"

# ╔═╡ c1a7bd4a-7a36-11eb-176a-f399eb6b5f49
# `CountingTropical{Float64}` has two fields `n` (number) and `c` (counting), both are of type `Float64`.
ein"ab,bc,cd,ad,de,ac->"([ising_bondtensor(CountingTropical{Float64}, J) for J in [J12, J23, J34, J14, J45, J13]]...)[]

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

# ╔═╡ 06bbead0-793f-11eb-0dec-c549b461b9cf
md"""
### Max 2-satisfiability problem
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
We define an energy function
```math
E = \sum\limits_{k=1}^{|C|} C_k(s_{i_k}, s_{j_k})
```
where ``C_k(s_{i_k}, s_{j_k})`` is ``+1`` if the ``k``th clause on $i$th and $j$th boolean variable is satisfied, otherwise, it is ``-1``. In the tensor network representation, we use a vertex tensor to represent a boolean variable

```math
T_{v}^n = \delta^{n,q=2}.
```
"""

# ╔═╡ 58f6f6eb-d722-4144-b091-5b6bd7f3e97c
function twosat_vertextensor(::Type{T}, n::Int) where T
	res = zeros(T, fill(2, n)...)
	res[1] = one(T)
	res[end] = one(T)
	return res
end

# ╔═╡ 528b4ef8-872a-11eb-14da-918dfd3110aa
md"If two vertices are related by a clause $C_k$, we add an bond tensor $T_k$ to connect two vertices.

```math
(T_{k})_{s_i s_j} = C_k(s_i,s_j)
```
"

# ╔═╡ fe9e0d86-ddab-4f39-a36f-5f42887780f6
function twosat_bondtensor(::Type{T}, src::Bool, dst::Bool) where T
	res = [T(1) T(1); T(1) T(1)]
	res[Int(src)+1, Int(dst)+1] = T(-1)
	return res
end

# ╔═╡ 268959f8-4162-4372-95b4-6f9fc47c23a6
md"The degeneracy can also be obtained easily with `einsum` (without using vertex tensors)"

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

# ╔═╡ 73f517de-7aed-11eb-03d1-db03dfb01a35
md"Since the resulting (counting) tropical number 11 is equal to the number of clauses, all clauses are satisfied, and the degeneracy is 16. (Note: this contraction is not computational efficient, for large scale problem, one should use the contraction order algorithms listed in Sec. 5.)"

# ╔═╡ 5f2243c4-793d-11eb-1add-392387bb559f
md"""
### Potts model

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
"""

# ╔═╡ d8daf729-c6fc-4d84-9cf4-bd4f3c6a3c15
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

# ╔═╡ c4c9caa6-872a-11eb-3c90-3b18e9762253
md"""and the vertex tensor

$$T_v^n=\delta^{n,q=3}$$
"""

# ╔═╡ 50c9910f-6b47-4480-9efd-768dae573b92
potts_vertextensor(T, q, n) = δtensor(T, q, n)

# ╔═╡ 344042b4-793d-11eb-3d6f-43eb2a4db9f4
md"""
### Maximum independent set
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
md"We can put a bond matrix at the bond to describe the independence restriction
```math
T_{e} = \begin{bmatrix}0 & 0 \\0 & -\infty\end{bmatrix}
```"

# ╔═╡ 64eb9dab-21bd-4412-8004-f82d7659ca2a
function mis_bondtensor(::Type{T}) where T
	res = ones(T, 2, 2)
	res[2, 2] = zero(T)
	return res
end

# ╔═╡ 06ebea54-872b-11eb-2903-157e0fc88de1
md"""
The vertex tensors are for counting the number of vertices
```math
(T_{v}^{n})_{s_i s_j} = \begin{cases}
 0, & s_i = s_j =\ldots s_n = 0\\
 1, & s_i = s_j =\ldots s_n = 1\\
 -\infty, &otherwise
\end{cases}
```
where $s_i,s_j,\ldots s_n \in \{0,1\}$.
"""

# ╔═╡ 0630d232-2791-4834-8076-3aba6c1deaee
function mis_vertextensor(::Type{T}, n::Int) where T
	res = zeros(T, fill(2, n)...)
	res[1] = one(T)
	res[end] = T(1)
	return res
end

# ╔═╡ 6800e29e-bee3-4670-b5b5-aaeee9e25046
md"We can build a tensor tensor network for solving this issue. The tensor network is contracted with a greedy order."

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
md"Or equivalently, we can transfer the weights on the vertex to bond tensors so that the vertex tensors can become δ tensors.
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

# ╔═╡ c2f987aa-7a36-11eb-0d03-4b6d328d8fa4
md"## Benchmarks

Note: Most of these benchmarks does not contract directly with `SimpleTensorNetworks` or `OMEinsum`. We mapped the tensor network to circuit simulations to achieve a better speed and smaller memory usage, see the last section or the following github repository.
"

# ╔═╡ d531c952-7ad9-11eb-1247-dd1913cc4678
html"""<div align="center"><a href="https://github.com/TensorBFS/TropicalTensors.jl">code is available on github <svg class="octicon octicon-mark-github v-align-middle" height="32" viewBox="0 0 16 16" version="1.1" width="32" aria-hidden="true"><path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"></path></svg></a></div>"""

# ╔═╡ 541f4062-7b22-11eb-2eb8-17585a3de9c3
md"**Ising spin glass on square lattices**"

# ╔═╡ 080bf23c-7ad8-11eb-37ac-01d2b6439f55
let
	img1 = viz_lattice(SquareLattice(32,32); r=0.01, ncolor=2)
	Compose.set_default_graphic_size(14cm*0.4, 7cm*0.8)
	leftright(updown(img1, md"We obtain the exact ground state energy of Ising spin glasses on square lattice up to $32^2$ spins."), updown(HTML("""<img src="https://user-images.githubusercontent.com/6257240/109566189-87bc5980-7ab1-11eb-9d08-99cd573007df.png" width=270px></img>"""), md"""Wall clock time for computing the ground state energy of the (a) Ising spin glass on an open square lattice with
``L^2`` spins. (tensor networks are contracted with [Yao.jl](https://github.com/QuantumBFS/Yao.jl))"""))
end

# ╔═╡ 64f18c2e-7b22-11eb-352f-9d6e228cef49
md"**Ising spin glass on cubic lattices**"

# ╔═╡ 9deb5b9a-7adf-11eb-3ba0-0d3716d7d603
let
	Compose.set_default_graphic_size(14cm*0.4, 7cm*0.8)
	θ = 2.31
	ϕ = 2.8
	cam_position = SVector(0.0, 0.0, 0.5)
	rot = RotY(θ)*RotX(ϕ)
	cam_transform = PerspectiveMap() ∘ inv(AffineMap(rot, rot*cam_position))
	Nx = Ny = Nz = 6
	nb = [nodestyle(:circle, fill(color); r=0.015) for color in ["#CC4040", "#4040CC"]]
	eb = bondstyle(:default, stroke("#333333"), linewidth(0.2mm); r=0.01)
	c = Cubic((0.05, 0.05, 0.05))
	x(i,j,k) = cam_transform(SVector(c[i-Nx/2-0.5,j-Ny/2-0.5,k-Nz/2-0.5])).data
	fig = canvas() do
		for i=1:Nx, j=1:Ny, k=1:Nz
			rand(nb) >> x(i,j,k)
			i!=Nx && eb >> (x(i,j,k), x(i+1,j,k))
			j!=Ny && eb >> (x(i,j,k), x(i,j+1,k))
			k!=Nz && eb >> (x(i,j,k), x(i,j,k+1))
		end
	end
	leftright(Compose.compose(context(0.5,0.5, 1.0, 1.0), fig), md"""Ground state and degeneracy on cubic lattice up to $6^3$ spins. 
""")
end

# ╔═╡ 72742d66-7b22-11eb-2fac-cb2534558248
md"**Ising spin glass on Chimera lattices**"

# ╔═╡ 4f3a12e0-7ad5-11eb-2b37-c95342185c3e
let
	Compose.set_default_graphic_size(6cm, 8cm)
	img = Compose.compose(context(0.0,0,8/6,1), viz_lattice(ChimeraLattice(8, 8); r=0.01, ncolor=2))
	leftright(updown(img, md"``\pm J`` Ising spin glass on the chimera graph of D-Wave quantum annealer of $512$ qubits in less than $100$ seconds and investigate the exact value of the residual entropy of $\pm J$ spin glasses on the chimera graph."), updown(html"""<img src="https://user-images.githubusercontent.com/6257240/109566350-bb977f00-7ab1-11eb-953f-127d7919e3e6.png" width=270px/>""", md"Wall clock time for computing the ground
state energy of Ising spin glass on the chimera graph with the
``L \times L`` unit cell (``8L^2`` spins). (tensor networks are contracted with [Yao.jl](https://github.com/QuantumBFS/Yao.jl))"))
end

# ╔═╡ e739e74c-7af0-11eb-104f-5f94da1bf0be
md"**3-state Potts model on square lattice**"

# ╔═╡ 80d6c2b6-7aef-11eb-1bf5-5d4f266dfa73
let
	Compose.set_default_graphic_size(5cm, 5cm)
	img = viz_lattice(SquareLattice(18,18); r=0.015, ncolor=3)
	leftright(updown(
		img,
		md"Ground-state energy, entropy, and computational
time of ``q = 3`` state Potts spin glass model on square lattices of sizes ``n = 4\sim 18``. Each data point is
averaged over ``100`` random instances computed on a single
GPU. As a comparison, the existing branch-and-cut method
with the Semi-Definite Programming energy lower bounds
method on the same model works up to ``9 \times 9`` lattices
(using 10 hours)"
		; width=300)
		,html"""<img src="https://user-images.githubusercontent.com/6257240/109578583-6d8c7680-7ac5-11eb-93eb-1b2748f2c90b.png" width=250px/>""")
end

# ╔═╡ 81e06ac6-7b22-11eb-3042-373a49bbdb49
md"**Ising spin glass and Max-2-SAT on random 3-regular graphs**"

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

	img = viz_graph(random_regular_graph(220, 3), linecolor="black", node_fontsize=0)
	#leftright(img, updown(html"""<img src="https://user-images.githubusercontent.com/6257240/109566350-bb977f00-7ab1-11eb-953f-127d7919e3e6.png" width=270px/>""", md"Wall clock time for computing the ground state energy of Ising spin glass on the chimera graph with the ``L \times L`` unit cell (``8L^2`` spins)."))
	Compose.set_default_graphic_size(7cm, 6cm)
	leftright(Compose.compose(context(0.0, 0.0, 6/7, 1.0), img), md"""
The spin glass on the random graphs: our method can compute optimal solutions and count the number of solutions for spin glasses and combinatorial optimization problems on on $3$ regular random graphs up to $220$ spins, on a single GPU. This is inaccessible by existing methods.
""")
end

# ╔═╡ 7bdf517e-79ff-11eb-38a3-49c02d94d943
md"## The Song Shan Lake Spring School (SSSS) Challenge"

# ╔═╡ 89d737b3-e72e-4d87-9ade-466a84491ac8
md"In 2019, Lei Wang, Pan Zhang, Roger and me released a challenge in the Song Shan Lake Spring School, the one gives the largest number of solutions to the challenge quiz can take a macbook home ([@LinuxDaFaHao](https://github.com/LinuxDaFaHao)). Students submitted more than 10 [solutions to the problem](https://github.com/QuantumBFS/SSSS/blob/master/Challenge.md). The quiz is"

# ╔═╡ 88e14ef2-7af1-11eb-23d6-b34b1eff8f87
md"""
In the Buckyball structure as shown in the figure, we attach an ising spin ``s_i=\pm 1`` on each vertex. The neighboring spins interact with an anti-ferromagnetic coupling of unit strength.

1. Get ``\ln Z/N``, where ``N`` is the number of vertices, and
```math
Z = \sum_{\{s\}}e^{-\sum_{i,j \in E}s_i s_j}
```

2. Count the ground state degeneracy.
"""

# ╔═╡ f867e24a-bbfa-44e7-819e-1051114f53f9
md"Now we are proud to announce a new solution to this quiz"

# ╔═╡ 5a5d4de6-7895-11eb-15c6-bda7a4342002
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

# ╔═╡ 9b1dc21a-7896-11eb-21f6-bfe9b4dc9ccf
let
	θ2 = 0.5
	ϕ2 = 0.8
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

# ╔═╡ b6560404-7b2d-11eb-21d7-a1e55609ebf7
# the positions of fullerene atoms
c60_xy = fullerene();

# ╔═╡ 6f649efc-7b2d-11eb-1e80-53d84ef98c13
# find edges: vertex pairs with square distance smaller than 5.
c60_edges = [(i=>j) for (i,(i2,j2,k2)) in enumerate(c60_xy), (j,(i1,j1,k1)) in enumerate(c60_xy) if i<j && (i2-i1)^2+(j2-j1)^2+(k2-k1)^2 < 5.0];

# ╔═╡ 20125640-79fd-11eb-1715-1d071cc6cf6c
md"The resulting tensor network contains 90 edge tensors and 60 vertex tensors."

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

# ╔═╡ 698a6dd0-7a0e-11eb-2766-1f0baa1317d2
md"Then we find a proper contraction order by greedy search"

# ╔═╡ 020cfb20-8228-11eb-2ee9-6de0fc7700b1
md"Seed for greedy search = $(@bind seed Slider(1:10000; show_value=true, default=42))"

# ╔═╡ ae92d828-7984-11eb-31c8-8b3f9a071c24
tcs, scs, c60_trees = (Random.seed!(seed); trees_greedy(c60_tnet; strategy="min_reduce"));

# ╔═╡ 2b899624-798c-11eb-20c4-fd5523f7abff
md"The resulting contraction order produces time complexity = $(round(log2sumexp2(tcs); sigdigits=4)), space complexity = $(round(maximum(scs); sigdigits=4))"

# ╔═╡ 8522456a-823c-11eb-3cc1-fb720f1cc470
SimpleTensorNetworks.contract(c60_tnet, c60_trees[]).array[]

# ╔═╡ 1c4b19d2-7b30-11eb-007b-ab03052b22d2
md"If you see a 16000 in the counting field, congratuations! The greedy contraction order can be visualized by dragging the slider (if you run it on your local host)"

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
	θ2 = 0.5
	ϕ2 = 0.8
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

# ╔═╡ 4c137484-7b30-11eb-2fb1-190d8beebbc3
md"""Since the complexity of tensor contraction is exponential to the number of legs involved, *"what is the optimal contraction order"* becomes one of the most important issues in tensor network contraction. The greedy algorithm we used here is efficient but not optimal. Finding the optimal contraction order itself is NP-hard.
"""

# ╔═╡ e302bd1c-7ab5-11eb-03f6-69dcbb817354
md"## Towards productivity

To solve large scale problems. You also need a good contraction order algorithm, e.g.
* Finding the optimal contracting order for `einsum`, [Classical Simulation of Quantum Supremacy Circuits](https://arxiv.org/abs/2005.06787).
* A heuristic algorithm for contraction order [Contracting Arbitrary Tensor Networks: General Approximate Algorithm and Applications in Graphical Models and Quantum Circuit Simulations](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.060503)

* For simple lattices, you can also map the problem to [Yao.jl](https://github.com/QuantumBFS/Yao.jl) for fast contraction.

The Tropical BLAS project is under the progress,

- TropicalGEMM on CPU for Tropical numbers ([nearly done](https://github.com/chriselrod/LoopVectorization.jl/issues/201)),
- TropicalGEMM on GPU (WIP),
- TropicalGEMM for CountingTropical numbers,
"

# ╔═╡ Cell order:
# ╟─c456b902-7959-11eb-03ba-dd14a2cd5758
# ╟─dfa8834c-e8c6-49b4-8bde-0816b573cbee
# ╟─121b4926-7aba-11eb-30e1-7b8edd4f0166
# ╠═5bb40ad6-7b33-11eb-0b31-63d5e47fa0e7
# ╟─3205a536-7a17-11eb-3473-b71305c96ca4
# ╟─592825a5-890a-4de8-984f-3d4ca48daca4
# ╟─5066553a-9b22-4f3b-b8ba-13a88291d9b9
# ╟─3208fd8a-7a17-11eb-35ce-4d6b141c1aff
# ╟─32116a92-7a17-11eb-228f-0713510d0348
# ╟─3217d1ca-7a17-11eb-14eb-a77ccfd983a0
# ╟─3221a326-7a17-11eb-0fe6-f75798a411b9
# ╟─32277c3a-7a17-11eb-3763-af68dbb81465
# ╟─322d2958-7a17-11eb-2deb-613b7680a5bb
# ╟─3237e33e-7a17-11eb-2869-b92d0801bc6e
# ╟─ec841be8-7a16-11eb-3337-376e26b7da25
# ╟─d0b54b76-7852-11eb-2398-0911380fa090
# ╠═af13e090-7852-11eb-21ae-8b94f25f1a4f
# ╠═d770f232-7864-11eb-0e9a-81528e359d39
# ╟─5d16a472-785b-11eb-1b94-dd6d8f860c24
# ╟─3372871c-785b-11eb-3092-4bbc419cb788
# ╠═2173a6cc-785b-11eb-1ab6-7fb875224dd9
# ╟─518b7d4e-785b-11eb-3b7c-1389065b9cbd
# ╠═2868b292-785b-11eb-015e-6b5613bd9e39
# ╟─98ae0960-797d-11eb-3646-c5b7e05d3f7c
# ╠═ca9a13b8-0272-4557-97d4-9168923d363f
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
# ╟─22118a36-7a36-11eb-18c3-dd2adac6118b
# ╠═c85217b6-7a23-11eb-04c5-fb4dc9c02ef1
# ╟─b3b986aa-7a1d-11eb-17d4-e5675015b221
# ╟─b5cd769e-7a1e-11eb-1d82-e1c265dfdd52
# ╟─6e1c507e-7a1a-11eb-05bc-dbccc3aebdf9
# ╟─b52ead96-7a2a-11eb-334f-e5e5ff5867e3
# ╠═624f57db-7f07-4281-a547-d229b9a8413a
# ╟─05109d30-7a29-11eb-320a-fb0b0d8e2632
# ╟─d0ecd3f2-7a2d-11eb-126d-7dab740d8e1f
# ╟─9c860e2a-7a2e-11eb-231f-63e9aca1daa0
# ╟─5efee244-7a2d-11eb-3782-b9d55086d623
# ╠═37472f2a-7a2a-11eb-1be3-13513d61fcb2
# ╟─023ebf7c-7b36-11eb-1c9f-430773395534
# ╟─96290770-7a20-11eb-0ac8-33a6492c7b12
# ╟─c1f90d6c-7a1d-11eb-2843-f971b5f6f3b0
# ╟─64e08a56-7a36-11eb-29fd-03662b4d6612
# ╠═b975680f-0b78-4178-861f-5da6d10327e4
# ╟─f54119ca-7a1e-11eb-1bec-bf855e34658d
# ╟─37102544-7abf-11eb-3ac4-6702dfc55425
# ╟─25553aa4-eef8-40d4-bcf6-0a27f50ab1da
# ╠═75764e5b-d2f2-4178-a00f-8f362b59a0b8
# ╠═15cf0c36-7a21-11eb-3e14-63950bcce943
# ╟─a4bf7de3-7794-4609-88e7-9092f496a7bb
# ╠═a9544e66-7a27-11eb-2b27-1d2124988fb2
# ╟─3bb2e0c2-7a28-11eb-1ea5-ab03d16bf0b3
# ╟─695e405c-786d-11eb-0a6e-bb776d9626ad
# ╟─1bb36c52-a171-4993-ac86-2250e1e87a01
# ╟─ade34905-5a61-4f71-b347-e02fab120b5d
# ╟─43101224-7ac5-11eb-104c-0323cf1813c5
# ╠═792df1aa-7a23-11eb-2991-196336246c43
# ╠═8388305c-7a23-11eb-1588-79c3c6ce9db9
# ╟─6faf5bf0-a8f5-4db0-839a-9389d541a690
# ╠═c1a7bd4a-7a36-11eb-176a-f399eb6b5f49
# ╟─4190393a-7ac4-11eb-3ac6-eb8e3574fdc9
# ╟─56fdb22c-7ac4-11eb-2831-a777d9ca89f3
# ╟─06bbead0-793f-11eb-0dec-c549b461b9cf
# ╟─ef2d2446-793f-11eb-223a-c5fe0ed5e367
# ╠═58f6f6eb-d722-4144-b091-5b6bd7f3e97c
# ╟─528b4ef8-872a-11eb-14da-918dfd3110aa
# ╠═fe9e0d86-ddab-4f39-a36f-5f42887780f6
# ╟─268959f8-4162-4372-95b4-6f9fc47c23a6
# ╠═b8c8999a-7aec-11eb-3ccd-69b48fcb93c2
# ╟─73f517de-7aed-11eb-03d1-db03dfb01a35
# ╟─5f2243c4-793d-11eb-1add-392387bb559f
# ╠═d8daf729-c6fc-4d84-9cf4-bd4f3c6a3c15
# ╟─c4c9caa6-872a-11eb-3c90-3b18e9762253
# ╠═50c9910f-6b47-4480-9efd-768dae573b92
# ╟─344042b4-793d-11eb-3d6f-43eb2a4db9f4
# ╟─80d764a8-7afd-11eb-3fb8-79169ca56c7e
# ╟─04d11828-7afa-11eb-3e73-1bbecf566f74
# ╟─5d95a598-7afa-11eb-10eb-db79fa44dd2a
# ╟─d29470d4-7afa-11eb-0afc-a34e39d49aa5
# ╟─75c37046-7b1b-11eb-00f5-7fc49f73f4d9
# ╠═64eb9dab-21bd-4412-8004-f82d7659ca2a
# ╟─06ebea54-872b-11eb-2903-157e0fc88de1
# ╠═0630d232-2791-4834-8076-3aba6c1deaee
# ╟─6800e29e-bee3-4670-b5b5-aaeee9e25046
# ╠═0405f4d8-7afb-11eb-2163-597b2edcf17e
# ╠═3a34aaec-7afb-11eb-1fc4-2fbc027753cc
# ╟─5f4e0fec-7afd-11eb-37c7-11b84027136a
# ╟─f13469bc-7afb-11eb-3dab-2b6cdf290f6f
# ╟─891c39b8-7b1b-11eb-2d70-67fd24021027
# ╟─5caa1b7c-7b21-11eb-394d-379351fe5170
# ╠═1b08b3ac-7b1e-11eb-2249-ddd787c549d4
# ╟─c2f987aa-7a36-11eb-0d03-4b6d328d8fa4
# ╟─d531c952-7ad9-11eb-1247-dd1913cc4678
# ╟─541f4062-7b22-11eb-2eb8-17585a3de9c3
# ╟─080bf23c-7ad8-11eb-37ac-01d2b6439f55
# ╟─64f18c2e-7b22-11eb-352f-9d6e228cef49
# ╟─9deb5b9a-7adf-11eb-3ba0-0d3716d7d603
# ╟─72742d66-7b22-11eb-2fac-cb2534558248
# ╟─4f3a12e0-7ad5-11eb-2b37-c95342185c3e
# ╟─e739e74c-7af0-11eb-104f-5f94da1bf0be
# ╟─80d6c2b6-7aef-11eb-1bf5-5d4f266dfa73
# ╟─81e06ac6-7b22-11eb-3042-373a49bbdb49
# ╟─e59d7a44-7ae7-11eb-3d93-3bc5cc46bc65
# ╟─7bdf517e-79ff-11eb-38a3-49c02d94d943
# ╟─89d737b3-e72e-4d87-9ade-466a84491ac8
# ╟─9b1dc21a-7896-11eb-21f6-bfe9b4dc9ccf
# ╟─88e14ef2-7af1-11eb-23d6-b34b1eff8f87
# ╟─f867e24a-bbfa-44e7-819e-1051114f53f9
# ╠═5a5d4de6-7895-11eb-15c6-bda7a4342002
# ╠═b6560404-7b2d-11eb-21d7-a1e55609ebf7
# ╠═6f649efc-7b2d-11eb-1e80-53d84ef98c13
# ╟─20125640-79fd-11eb-1715-1d071cc6cf6c
# ╠═c26b5bb6-7984-11eb-18fe-2b6a524f5c85
# ╟─698a6dd0-7a0e-11eb-2766-1f0baa1317d2
# ╟─020cfb20-8228-11eb-2ee9-6de0fc7700b1
# ╠═ae92d828-7984-11eb-31c8-8b3f9a071c24
# ╟─2b899624-798c-11eb-20c4-fd5523f7abff
# ╠═3674a622-823b-11eb-3991-d5771010237b
# ╠═94b870d2-8235-11eb-33e7-35bf5132efd6
# ╠═8522456a-823c-11eb-3cc1-fb720f1cc470
# ╟─1c4b19d2-7b30-11eb-007b-ab03052b22d2
# ╟─58e38656-7b2e-11eb-3c70-25a919f9926a
# ╟─12740186-7b2f-11eb-35e4-01e6f9ffbb4d
# ╟─c1c74e70-7b2c-11eb-2f26-21f54ad00fb2
# ╟─4c137484-7b30-11eb-2fb1-190d8beebbc3
# ╟─e302bd1c-7ab5-11eb-03f6-69dcbb817354
