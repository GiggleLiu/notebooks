### A Pluto.jl notebook ###
# v0.17.2

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 4a96f5c9-37b4-4a8a-a6bd-8a4b4440eb49
using Pkg; Pkg.activate()

# ╔═╡ 2a145cba-26b0-43bd-9ab2-13818d246eae
using Revise, PlutoUI, Viznet, Compose

# ╔═╡ a8b06352-5265-4b89-bd82-b31f3cdac391
using Latexify

# ╔═╡ 57a3616e-49af-40b7-a000-a4ecc81af84e
using BitBasis

# ╔═╡ 40563837-b1e9-4df9-9cc7-8841d0068973
using Yao, YaoPlots

# ╔═╡ 675b3398-01d7-4949-bb0e-7cdf9b805c69
using SymEngine: Basic  # import the symbolic data type

# ╔═╡ 0f8d63ac-f677-4889-b033-2a93f62be700
using YaoExtensions: variational_circuit

# ╔═╡ 342685aa-5159-11ec-13fd-fb8954106bca
using Plots

# ╔═╡ 5fbddde1-3594-43c4-9f56-c0ae258d926f
using Yao.ConstGate: P0

# ╔═╡ 6f131f36-0b4e-4570-8527-620297fae48e
using YaoToEinsum, OMEinsumContractionOrders, OMEinsum

# ╔═╡ e11c26c0-e534-45fc-bb1c-c0f2ce4016db
SPACE = html"&nbsp; &nbsp; &nbsp; &nbsp;"

# ╔═╡ e6306a69-bd6a-4c01-9c6e-1cb523668019
begin
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
	
	function highlight(str)
		HTML("""<span style="background-color:yellow">$(str)</span>""")
	end
end;

# ╔═╡ 86b15fb0-f112-4689-8663-8cc6c0a8fb2a
html"<button onclick='present();'>present</button>"

# ╔═╡ 0f099c85-f039-477e-a70d-a3801cbb2656
let
	img = html"""<img src="https://swarma.org/wp-content/uploads/2020/01/2020011414445650.png" width=200/>"""
	md"""
# 量子线路与量子过程的经典模拟
主讲人： 刘金国
### Goals
*  What is quantum computing and $(highlight("Full amplitude")) and $(highlight("tensor network")) based quantum circuit simulation.
* Introduce some $(highlight("open source packages")) for quantum simulation in Julia

$img
	
# Contents
* From classical adder to reversible adder,
* What is new in Quantum, the Deutsch-Jozsa algorithm,
* Towards faster simulation, a tensor network based quantum circuit simulation (basic).

"""
end

# ╔═╡ 48854a73-4896-4542-9ad4-15ae87418f1d
md"# Background knowledge"

# ╔═╡ c1d40103-1710-4221-b414-0958c13fb95f
md"""
I wish you are familiar with the following notations (~ 3rd year undergraduate)
* ``|\psi\rangle`` is a quantum state, or "ket",
* ``H`` is a Hamiltonian, it determines the dynamics of the quantum state as ``|\psi(t)\rangle = e^{-iHt}|\psi(0)\rangle``,
* ``\mathcal{O}`` is an observable, known as a Hermitian operator. One can measure this observable. The expectation value is the observed result is $\langle\mathcal{O}\rangle = \langle\psi|\mathcal{O}|\psi\rangle$
* ``X, Y`` and ``Z`` are Pauli operators.
"""

# ╔═╡ 4c8a4cea-63f4-49f0-82b9-dfa608be46bf
md"## If not..."

# ╔═╡ eb75810a-a746-4e4a-889f-c87c8d1e153f
html"""
<p>
<img src="https://static.docsity.com/documents_first_pages/notas/2012/04/17/dfe728fbf6999820eecb2c8fc0e773b6.png" width=300 style="vertical-align:sub;"><img src="https://images-na.ssl-images-amazon.com/images/I/51X+dIBIeZL._SY344_BO1,204,203,200_.jpg" width=290 style="vertical-align:super;">
</p>
"""

# ╔═╡ 330356f0-d7d7-425f-b16c-d1e9f97494b5
md"
* Physics background: Sakurai -> Neilsen
* Computer science background: Sakurai <- Neilsen
"

# ╔═╡ b83f7675-68fa-44b5-8681-c85984eeb877
md"# From classical computation to reversible computation"

# ╔═╡ b7d90c25-2f66-4ace-8d52-841ea376b3f9
md"## How classical gates work"

# ╔═╡ 2e9c5b7f-14dc-4ea1-a286-6a68c3e96b3e
md"Truth table of Logic gates"

# ╔═╡ 59c24944-9a8d-4a7f-8b0b-c08a000cd655
md"""
![](https://welkersclass.weebly.com/uploads/8/5/5/8/85585096/2479870_1_orig.gif)
"""

# ╔═╡ 2c6d024a-3187-4976-af27-393af8826a2d
md"This is an adder"

# ╔═╡ 165b369a-3d0a-4ae7-99f1-de0297f93707
md"![](http://www.worldofindie.co.uk/wp-content/uploads/2018/02/4-bit-adder-subtractor-circuit.png)"

# ╔═╡ c841f2e2-b907-4f74-be57-968ca339bec4
md"
FULL_ADDER(cin, a, b) -> cout, c
```
000 → 00
001 → 01
010 → 01
011 → 10
100 → 01
101 → 10
110 → 10
111 → 11
```
"

# ╔═╡ e1b6b7a9-2d26-4f43-a0bd-54f7ad22a5b3
md"""## $(highlight("Univeral")) gate or gate set"""

# ╔═╡ 16febe07-7f52-4b74-b9d6-e20fd8b05ab3
md"They are univeral:

* \{NAND\}
* \{NOR\}
* \{NOT, AND\}
* \{NOT, OR\}
"

# ╔═╡ d631887a-a222-402d-94bb-ecc6db6bea56
md"### Example"

# ╔═╡ 42f3e654-1836-491a-8119-b03b93822f45
md"NAND = 
```
00 → 1
01 → 1
10 → 1
11 → 0
```
"

# ╔═╡ cc0d4b1c-07b6-42ed-b1e2-80581a02ee6f
md"NOT(x) = NAND(x, x)"

# ╔═╡ c5694b66-d023-42c1-ae62-b7218ba8ebe7
md"AND(x ,y) = NOT(NAND(x, y))"

# ╔═╡ c5812ca6-c4ca-4211-9d2d-df498fd7a2da
md"OR(x, y) = NAND(NOT(x), NOT(y))"

# ╔═╡ 6b4e2611-5a5a-4646-aa8b-fec0d84d240c
md"# Reversible Gates"

# ╔═╡ ad47b86e-d8da-401d-a1a2-1139664adee5
md"""
1. number of inputs equals the number of outputs
2. the truth table is permutation
"""

# ╔═╡ ed0dc5f3-0da3-4d4f-9e3b-d4f471aa2e03
md"![](https://www.mathsisfun.com/sets/images/function-mapping.svg)"

# ╔═╡ ec1613ad-482e-4c68-b130-b51725f8d94e
md"""A reversible computer is $(highlight("equally powerful")) as an irreversible one"""

# ╔═╡ bde5a21b-7db3-40c2-8c00-54182987dcfd
md"## NOT gate -> X gate"

# ╔═╡ 07f36b6e-a98d-4bee-9091-b1a717dab0e8
not_gate(n, i) = put(n, i=>X);

# ╔═╡ a4bd7a89-7dda-4652-a30d-3e2b02079b1e
vizcircuit(not_gate(1, 1))

# ╔═╡ 91f767ea-5afa-4e04-9101-d1c9b45a5b1e
md"""
```
0 → 1
1 → 0
```
"""

# ╔═╡ e188093d-1ef1-4e70-966b-77cc0761a801
vmat(prefix, gate) = Latexify.LaTeXString(prefix * latexify(mat(Basic, gate); env=:raw).s);

# ╔═╡ f6533ed1-b419-43a7-b03e-3a538d2e46da
md"## A compact way of representing truth table"

# ╔═╡ cac05229-df74-4679-89ab-e65003b7d773
md"""
A state is a one hot vector, state 0:
```math
|0\rangle = \left[\begin{matrix}1 \\ 0\end{matrix}\right]
```

state 1:
```math
|1\rangle = \left[\begin{matrix}0 \\ 1\end{matrix}\right]
```
"""

# ╔═╡ 4518fcfa-9f04-4c2c-82b8-719033f20ac9
md"""

A **classical** **reversible** gate is an permutation matrix multiplied on the input state.
"""

# ╔═╡ bded6dc0-bce0-404e-a1f2-36712f26f2d6
vmat("X = ", X)

# ╔═╡ beed43ce-b24c-4edd-9da2-30bdd72c9411
md"""
```math
|\psi\rangle_{out} = X |\psi\rangle_{in}
```
"""

# ╔═╡ 0c267c28-967b-4268-8372-00d7b48b7b8e
md"## XOR gate -> Controlled not gate, or CNOT gate"

# ╔═╡ 7440555d-685f-4447-a44a-0463b37bd43c
xor_gate(n, i, j) = control(n, i, j=>X);

# ╔═╡ fed37413-c341-4dcf-a637-404d2f186b2b
vizcircuit(xor_gate(2, 1, 2); show_ending_bar=false)

# ╔═╡ 69d1fc34-3e02-4c4b-83a7-6d05a0d69405
md"""
```
00 → 00
01 → 01
10 → 11
11 → 10
```
"""

# ╔═╡ 67ad78c7-930a-48b1-b8d3-d9e549ee7379
vmat("CNOT = ", xor_gate(2,1,2))

# ╔═╡ 523fdddf-e76c-46f0-bda0-4d3e19c7a82d
md"""input state is a $(highlight("onehot vector"))
```math
|00\rangle = \left[\begin{matrix}1\\0\\0\\0\end{matrix}\right]
```
```math
|01\rangle = \left[\begin{matrix}0\\1\\0\\0\end{matrix}\right]
```
```math
|10\rangle = \left[\begin{matrix}0\\0\\1\\0\end{matrix}\right]
```
```math
|11\rangle = \left[\begin{matrix}0\\0\\0\\1\end{matrix}\right]
```
"""

# ╔═╡ f311befa-ba09-42e6-ac9b-b59450162ebd
md"## AND gate -> Toffoli gate"

# ╔═╡ 816ebc7b-43ec-47e2-a24a-86ac70dd6afe
and_gate(n, i, j, k) = control(n, (i, j), k=>X);

# ╔═╡ 893baeda-d5b0-43c1-9018-8d18c06486ca
vizcircuit(and_gate(3, 1, 2, 3); show_ending_bar=false)

# ╔═╡ 9a5227f6-0db8-43c1-9b79-80897faa86d0
md"""
```
000 → 000
001 → 001
010 → 010
011 → 011
100 → 100
101 → 101
110 → 111
111 → 110
```
"""

# ╔═╡ cee8db91-7e84-4728-aef6-ea861c62ff96
vizcircuit(control(3, (1,2), 3=>X); show_ending_bar=false, starting_texts=[1,1,0], starting_offset=-0.3, ending_texts=[1,1,1], ending_offset=0.3)

# ╔═╡ 38194214-3bf1-4229-9fe8-37282b30a5ad
vmat("Toffoli = ", control(3, (3,2), 1=>X))

# ╔═╡ 4591b3a0-ff3f-46c8-8c8d-ef7c6ba76cc5
md"## NAND gate"

# ╔═╡ 7c1cd22b-c478-4dcf-9828-ae6aad3df95a
nand_gate(n, i, j, k) = chain(control(n, (i, j), k=>X), put(n, k=>X));

# ╔═╡ ae12828c-350a-4251-96cf-5ac2dc0fc0ee
vizcircuit(nand_gate(3,1,2,3))

# ╔═╡ 6cc597c9-e370-48ea-a5c5-422904d4d8a0
md"## Or gate"

# ╔═╡ b2bf15bb-1350-4304-a932-c87b09558115
or_gate(n, i, j, k) = chain(kron(n, i=>X, j=>X, k=>X), control(n, (i, j), k=>X), kron(n, i=>X, j=>X));

# ╔═╡ 0e63c1eb-fea6-411a-bae6-0bc90dee6bc7
vizcircuit(or_gate(3, 1, 2, 3); show_ending_bar=false)

# ╔═╡ 0b42d7bb-b9cf-42ef-aad7-5fc5d5918be3
md"## Reversible Full Adder"

# ╔═╡ 622fdaac-87dd-4c22-85b9-470510566480
# `s`, `cout` and ancillas (`x` and `y`) are initialized to 0
# flush data in `cin`
full_adder(n, a, b, cin, cout, x, y) = chain(and_gate(n, cin, a, x), xor_gate(n, a, cin), and_gate(n, b, cin, y), or_gate(n, x, y, cout), xor_gate(n, b, cin));

# ╔═╡ 26e5d619-7119-49bc-8907-17ae0db424f5
vizcircuit(full_adder(6, 1:6...); scale=0.7, starting_texts=["a", "b", "cᵢ", "0", "0", "0"], ending_texts=["a", "b", "cᵢ'", "cₒ", "?", "?"])

# ╔═╡ a7d85b82-a705-4f3b-a371-06a87071335d
function add_circuit(n::Int)
	nbit = 5n+1
	as = 1:n
	bs = n+1:2n
	cin = 2n+1
	carries = 2n+2:3n+1
	xs = 3n+2:4n+1
	ys = 4n+2:nbit
	c = chain(nbit)
	cs = zeros(Int, n+1)
	for i=1:n
		cout = carries[i]
		blk = full_adder(nbit, as[i], bs[i], cin, cout, xs[i], ys[i])
		push!(c, blk)
		cs[i] = cin
		cin = cout
	end
	cs[end] = cin
	return c, cs
end;

# ╔═╡ 86b398d0-844d-4d06-9ede-79d618502cce
let
	n = 2
	c, outputs = add_circuit(n)
	N = nqubits(c)
	m = Yao.Measure(N; locs=outputs)
	mc = chain(c, m)
	#vizcircuit(mc; scale=0.5)
	product_state(bit"00000001111") |> mc
	m.results
end;

# ╔═╡ e2e4ec60-2737-4560-89b1-1e14a35044e8
function calculate_binaryadd(n::Int, x, y)
	c, outputs = add_circuit(n)
	res = product_state(bit_literal([takebit(x, i) for i=1:n]..., [takebit(y, i) for i=1:n]..., zeros(Int, 3*n+1)...)) |> c |> measure!
	Int(readbit(res, outputs...))
end;

# ╔═╡ aea490af-0bdd-4930-9ad2-7d9a13e08c46
calculate_binaryadd(2, 2, 2);

# ╔═╡ 7f84eed8-edd8-4f8e-a2b1-3ad862285934
md"## A 4 bit binary adder"

# ╔═╡ c0078012-3d81-4584-a050-9a58802d08a9
gatecount(add_circuit(4)[1]);

# ╔═╡ ea458fa2-1f9f-46e1-88da-942034d0fa73
let
	n = 4
	circuit, out = add_circuit(n)
	vizcircuit(circuit; scale=0.3, starting_texts=[["a$i" for i=1:n]..., ["b$i" for i=1:n]..., ["0" for i=1:3n+1]...], ending_texts=[["a$i" for i=1:n]..., ["b$i" for i=1:n]..., [i+2n ∈ out ? "c$(i)" : "?" for i=1:3n+1]...])
end

# ╔═╡ 3381001d-1120-4b88-ac01-5ca861f0a9be
md"Note: ancillas must be initialized to state 0"

# ╔═╡ 7524029c-2b7e-465e-9827-c993d6cdd34a
Yao.YaoBlocks.LuxurySparse.SparseMatrixCSC(mat(add_circuit(4)[1]))

# ╔═╡ 22cd8a36-c76f-4a8c-a8a2-e1924136012a
md"Adder with uncomputing"

# ╔═╡ 051548ff-dd39-4d55-ac53-e8e2bacec68e
function adder_with_uncomputing(n)
	add, outputs = add_circuit(n)
	n1 = nqubits(add)
	chain(subroutine(n1+n, add, 1:n1), [control(n1+n, outputs[i], n1+i=>X) for i=1:n]..., subroutine(n1+n, add', 1:n1)), collect(n1+1:n1+n)
end;

# ╔═╡ e7b2e327-952e-4062-9193-19653eaee19c
let
	n = 4
	circuit, out = adder_with_uncomputing(n)
	vizcircuit(circuit; scale=0.3, starting_texts=[["a$i" for i=1:n]..., ["b$i" for i=1:n]..., ["0" for i=1:4n+1]...], ending_texts=[["a$i" for i=1:n]..., ["b$i" for i=1:n]..., [i+2n ∈ out ? "c$(i)" : "0" for i=1:4n+1]...])
end

# ╔═╡ cb2c496c-e9dc-4666-b41b-d97cca377047
md"## Universal reversible gate

* {Toffoli}
* {Fredkin}
"

# ╔═╡ 78cb38db-b0b6-49f7-bcaf-1a2df58b8150
md"CNOT = "

# ╔═╡ 334db854-7805-4548-a414-ad7d215fe387
vizcircuit(and_gate(3, 1, 2, 3); starting_texts=["1", "a", "b"])

# ╔═╡ 40956c0e-3833-4397-a4a8-c46e1890ef39
md"NOT = "

# ╔═╡ 8bebf527-d197-48b0-9778-f71b0a5dad77
vizcircuit(and_gate(3, 1, 2, 3); starting_texts=["1", "1", "a"])

# ╔═╡ ca094f81-bef3-40eb-9279-086fe7eb506a
md"## Irreversible v.s. Reversible
* Reversible computing can be much more energy efficient (arxiv: 1803.02789), because erasing 1 bit information costs a least ``k_bT`` energy (known as the Landauer's principle).
* Reversible computing may have polynomial time/space overhead.
"

# ╔═╡ 0840e7ea-aa63-44df-a788-ad18ac842006
md"## Balanced or constant?"

# ╔═╡ e5004508-2b9a-4acd-91ba-ad21fc7d1b16
md"""
``f: \{0, 1\}^n \rightarrow \{0, 1\}``

* *Balanced*: $f$ has an equal probability to output vaue 0 and 1.
* *Constant*: $f$ is a constant function.
"""

# ╔═╡ 5c03fc94-ce6d-490b-a3ee-911692877e6a
md"Requires: You can not execute the function twice!"

# ╔═╡ 2da83f06-e639-4f75-a162-cb561a8207ca
constantf(n) = chain(control(n+1, 1, n+1=>X), put(n+1, n+1=>X), control(n+1, 1, n+1=>X));

# ╔═╡ 960b0cdc-a769-4f03-8842-4e9d86e67cef
md"e.g. A constant function"

# ╔═╡ 408f534d-5553-4a95-8d14-cabf36f8aff5
vizcircuit(constantf(4), starting_texts=[["a$i" for i=1:4]..., "0"], ending_texts=[["a$i" for i=1:4]..., "f"])

# ╔═╡ 41244f02-796b-4afb-9825-0b21d216b78b
md"e.g. A balanced function
```
f(x, y) = (x + y) % 2
```
"

# ╔═╡ b8069c02-d422-4162-ad9c-c4b94870db6d
md"""# Summary 1
* Classical Adder -> Reversible Adder,
* Reversible circuit notation,
    * NOT gate
    * AND gate
    * XOR gate
* One hot representation of a classical state, which has a probabilistic interpretation. A reversible gate is a permutation matrix.
* Univeral gate sets
* The problem of distinguishing balanced or constant functions.

### TODO
* A tutorial of Yao,
* Build a quantum adder with Yao,
* Use a quantum algorithm to solve the "balanced or constant" problem,
* Tensor network based simulation of quantum circuits,
"""

# ╔═╡ 0f655562-0c18-43a8-83f0-4d11b7eeb068
vizcircuit_with_inputs(g, inputs) = vizcircuit(g; show_ending_bar=false, starting_texts=inputs, starting_offset=-0.3, ending_texts=measure!(product_state(bit_literal(inputs...)) |> g), ending_offset=0.3);

# ╔═╡ 774def48-8f9f-425a-aa47-ffe396692247
md"# Quantum bits and quantum gates"

# ╔═╡ 8017edf3-05b2-4cee-8a4e-90331a681037
let
	img = html"""<img src=https://docs.julialang.org/en/v1/assets/logo.svg width=100>"""
	md"""
## Why Julia?
$(img) is a **high-level**, **high-performance**, **dynamic programming language**. While it is a general-purpose language and can be used to write any application, many of its features are well suited for **numerical analysis** and **computational science**.

![](https://6lli539m39y3hpkelqsm3c2fg-wpengine.netdna-ssl.com/wp-content/uploads/2021/01/Julia-Adoption.png)
"""
end

# ╔═╡ d11e3541-6e03-41d7-a04b-fcce244058f7
md"## Julia 中文社区2021 meetup 会议视频 （Bilibili）"

# ╔═╡ 33bd9c30-8570-4580-94a7-adbfec13bedc
html"""<img src="https://user-images.githubusercontent.com/6257240/146650560-c6085dd3-03ea-4681-9f7f-6d9d80853922.png" width=500/>"""

# ╔═╡ 40a3bc85-be6d-4d1b-8376-da5a33c7057b
md"感谢陈久宁以及其它JuliaCN核心成员的努力以及集智的宣传！"

# ╔═╡ 52d14933-99d4-4ee6-9eff-7be4a5722334
md"""## Why Yao?

$(html"<img src=https://camo.githubusercontent.com/477280de44a6d4408d3a3255d3d82a615a27ac1c5120063ef7d6b2f6640befb8/68747470733a2f2f79616f7175616e74756d2e6f72672f6173736574732f6c6f676f2e706e67 width=100>")
[https://yaoquantum.org/](https://yaoquantum.org/) (arXiv:1912.10877)

Yao is an open source framework that aims to empower quantum information research with software tools. It is designed with following in mind:

* quantum algorithm design;
* quantum software 2.0;
* $(highlight("quantum computation education."))

by Roger Luo and Jinguo Liu, under the guidance of Lei Wang and Pan Zhang.

The most popular quantum simulator in Julia,
one of **Top 50** most popular Julia packages. It provides

* State of the art performance, with GPU support,
![](https://github.com/Roger-luo/quantum-benchmarks/raw/master/images/pcircuit.png)
* Matrix representation for operators,
* Built-in automatic differentiation engine,
* Generic data type (symbolic computation, tropical algera),
"""

# ╔═╡ e5d6f98c-aa84-48cf-bbca-cd72333d67f7
md"## Quantum state and Born's rule"

# ╔═╡ 99b9d334-40ed-4636-b299-1afe40402a12
md"Onehot vector -> Quantum state"

# ╔═╡ 400f697d-fbdb-4067-881f-e004f50bbc0f
md"""A quantum state is represented as a $(highlight("normalized complex valued")) vector. e.g. the following vector represents a quantum system with 4 states (``s=1,2,\ldots,4``)"""

# ╔═╡ b3daf0a3-4d1e-4474-bb99-948c718059ce
reg_random = Yao.rand_state(3)

# ╔═╡ 30c16e7d-417c-4098-a879-0914d711ff56
Yao.statevec(reg_random)

# ╔═╡ 6a719cf5-d8a4-4369-95e7-c1dfbc5cebf8
md"""
Its elements ``\langle s|\psi\rangle =\sqrt{p_s}e^{i \phi_s}``. Where $p$ is the probability given by Born's rule.
"""

# ╔═╡ 6f64057b-130d-44c1-8957-e4c7c56f7cb4
reg_product = Yao.product_state(bit"010")

# ╔═╡ 417b6e83-4320-40eb-a91c-9e6cf5f36be8
md"product state corresponds to the one hot vector."

# ╔═╡ cd39c31a-6fc9-489f-8045-cd3eca229a76
Yao.statevec(reg_product)

# ╔═╡ 8d1c023b-4ee6-4201-99f6-0245cee70cda
reg_zero = Yao.zero_state(3)

# ╔═╡ 7278c646-3c7b-4b18-a989-cb98abefd5cc
Yao.statevec(reg_zero)

# ╔═╡ 30100576-9fa2-4825-bfa6-35ae6cb920bb
md"# Quantum Gates and Operators"

# ╔═╡ 217e42d4-d226-42b2-bbd7-de85bc467ec4
md"## Quantum gate is time evolution operator"

# ╔═╡ 0e1e5385-db28-474f-a14a-cf3fd1efcd94
md"permutation matrix -> unitary matrix"

# ╔═╡ 681cdc3b-348f-45c4-9031-6298e81d995d
md"""
Hamiltonian is a **Hermitian** matrix
"""

# ╔═╡ b3813f2d-4e9d-4789-a653-d28fd1b98eda
hami = kron(X, X) + kron(Y, Y) + kron(Z, Z)

# ╔═╡ a9d19c74-1b83-4211-9c83-95ad1059d432
mat(hami)   # get the matrix representation

# ╔═╡ 48427378-face-48eb-b635-f3f460a51e26
ishermitian(hami)

# ╔═╡ 45d19896-ae93-4f9a-b1fe-bb143c7ff580
md"Time evolve
```math
|\psi(t)\rangle = e^{-iHt} |\psi(0)\rangle
```
"

# ╔═╡ bf983162-66e3-40e0-afce-7f5b7728c490
let
	gadget = @bind evolve_time Slider(0:0.03:2; show_value=true)
	md"evolve_time = $gadget"
end

# ╔═╡ bf357f44-54e2-49c0-b0b6-263abb1dbcb0
gate = time_evolve(hami, evolve_time)

# ╔═╡ 913bc173-8cfd-4680-ad13-36ea806acecb
mat(gate)

# ╔═╡ 09df3c8e-3d25-4ba2-bf19-9947224dfca4
mat(Basic, SWAP)

# ╔═╡ 45ace908-90df-488c-9258-459c0344baae
isunitary(gate)

# ╔═╡ de21a067-4587-461d-94bb-a295b6768f96
ψ0 = rand_state(2)

# ╔═╡ d85a2028-9b6f-4d6c-ba61-fecef3fd4876
ψt = apply(ψ0, gate)  # applying gate on state ψ0

# ╔═╡ 4a7e56ac-a69f-4b1f-b585-3a92f43e42e8
statevec(ψt)

# ╔═╡ f4f18b91-86f6-460c-b5de-4c52931a4098
md"## Quantum universal gate set
* \{Toffoli, H\}

* \{CNOT, H, S, T\}
"

# ╔═╡ bfb9fcb5-2555-490a-9dcc-48cc91f8ab3c
md"""
## Primitive Operators
"""

# ╔═╡ 2e17662d-24ed-40c9-93b0-4cef526c3a75
@bind selected_gate Select([X, Y, Z, I2, ConstGate.P0, ConstGate.P1, ConstGate.Pu, ConstGate.Pd, ConstGate.T, ConstGate.S, SWAP, ConstGate.CNOT, ConstGate.Toffoli, Yao.Measure(1)])

# ╔═╡ 21b11fac-5efb-4fd5-a1b0-e684d215a46c
vizcircuit(selected_gate)

# ╔═╡ 3fa73e8a-fe46-475f-8346-e54f52c144f4
md"Hermitian: ``\mathcal{O} = \mathcal{O}^\dagger``"

# ╔═╡ 31c1c4fc-7476-4066-8462-bf18f8d69966
ishermitian(selected_gate)

# ╔═╡ 81a55867-4158-4371-a12c-b08e3d64c0cb
md"Unitary: ``\mathcal{O}^\dagger\mathcal{O} = I``"

# ╔═╡ 7f15cbb0-09be-46ee-852c-e43dde9bc4f5
isunitary(selected_gate)

# ╔═╡ cc4e7d7b-0596-4bc0-b23e-f6861fcd5260
md"Reflexive: ``\mathcal{O}^2 = I``"

# ╔═╡ 8b03f2f1-f38c-4d5e-b778-0c3b4aaf910d
isreflexive(selected_gate)

# ╔═╡ a2129e65-7e73-4b42-9924-e88d60893ed2
Matrix(mat(Basic, selected_gate))

# ╔═╡ 4dd7f18c-2c0f-4405-bbfe-706699ea8958
md"## Rotation gate"

# ╔═╡ faf2b47d-ab3b-4087-8cda-a564303c3bc9
md"``R_g(\theta) = e^{-ig\theta/2}``, `g` is a reflexive and hermitian operator."

# ╔═╡ 8b1606b5-78d8-4020-8597-1cf9f2364a2b
mat(rot(X, Basic(:θ)))

# ╔═╡ 04693ee8-d0b2-4ffa-80b2-40b9149eceea
vizcircuit(rot(X, Basic(:θ)))

# ╔═╡ 8be8d92a-d882-4e77-a3c3-9846babeda2b
mat(rot(SWAP, Basic(:θ)))

# ╔═╡ 597bff8b-f7fa-4fc8-bc59-f06c934f211e
md"""
## Put a gate at a qubit
"""

# ╔═╡ 421abcef-a1be-4297-a159-adcd92332d6b
md"
`put(n, (i,j...,)=>G)`
```math
I_1\otimes I_2 \otimes \ldots G_{i,j,\ldots} \ldots \otimes I_n
```
"

# ╔═╡ 3388350c-fd8b-4cae-bf34-71446bd747ca
g_put1 = put(2, 2=>X)

# ╔═╡ 99847765-fd49-4c9d-914b-5d47eac89490
vizcircuit(g_put1)

# ╔═╡ fbd1ea6e-16da-4289-b0e3-d8e8a5a55541
mat(Basic, g_put1)

# ╔═╡ 5dafbbbd-d6f7-423c-8a46-41b1596f83a3
g_put2 = put(3, (3, 1)=>SWAP)

# ╔═╡ 40889c09-1a10-4829-9a4e-5dc46d9045ef
vizcircuit(g_put2)

# ╔═╡ 984db7f2-30e0-421d-9f27-874553b6fae2
mat(Basic, g_put2)

# ╔═╡ 44d6b7ba-43e6-43f7-bbed-a727edd22e0c
md"## Control gates"

# ╔═╡ f35b9789-788b-450d-af8e-472412efce77
md"`control(n, (c,d,...), (i,j,...)=>G)`
```math
\begin{align}
&I_1\otimes I_2\otimes \ldots P_{c,d,\ldots} \ldots G_{i,j,\ldots} \ldots \otimes I_n + I_1\otimes I_2\otimes \ldots (1-P_{c,d,\ldots}) \ldots I_{i,j,\ldots} \ldots \otimes I_n\\
&P_{c,d,\ldots} = P_c \otimes P_d \otimes \ldots\\
&P_{c} = \begin{cases}
|1\rangle\langle 1|_{|c|}, c > 0\\
|0\rangle\langle 0|_{|c|}, \text{otherwise}
\end{cases}
\end{align}
```"

# ╔═╡ 34a0b380-ac20-4d88-b498-969745a0df7c
g_ctrl1 = control(2, 2, 1=>X)

# ╔═╡ 3aaa7e49-9507-4384-ac7e-f4d486127811
vizcircuit(g_ctrl1)

# ╔═╡ ac2f9c18-f4bc-4077-9ad1-1f4f0d75a7b1
md"Inverse control"

# ╔═╡ 47da74c7-0f24-4dd5-b6f7-88c22a81de5a
g_ctrl2 = control(2, -2, 1=>X)

# ╔═╡ e06b23aa-7fae-4495-9b69-43cea3682a28
vizcircuit(g_ctrl2)

# ╔═╡ 915185d6-3b36-437e-bac0-22cd5503bba2
vizcircuit(control(3, (3,2), 1=>X))

# ╔═╡ 67d9d529-aeec-4ede-a0f7-85dc7f4fba5a
md"## Compose multiple gates"

# ╔═╡ 0a6d5a41-f42b-42e7-adf7-93a40f67c9f3
md"`chain(G1, G2, ..., Gn)` or `Gn * ... * G_2 * G_1`
```math
G_n G_{n-1} \ldots G_1
```
"

# ╔═╡ dc46cbbe-0eb5-41f9-ad5f-f8a7fba5ff9c
mat(Basic, chain(X, Y))

# ╔═╡ 0a7cab72-a50b-4d58-877a-4a9d939b198f
g_chain2 = chain(put(2, 1=>X), control(2, 1, 2=>X))

# ╔═╡ 6b2d8bb9-5a4f-4699-ac92-bc8a4003918f
vizcircuit(g_chain2)

# ╔═╡ 17a574c9-664a-4cba-8571-60a2551f37f7
md"## Add up terms to obtain Hamiltonians"

# ╔═╡ 3a89d593-ae7e-454c-bb18-0baa4cc06a25
md"+(H1, H2, ..., Hn)
```math
H_n + H_{n-1} \ldots H_1
```
"

# ╔═╡ 4559f70b-5a9b-478d-ba9a-3a35f6fc9716
+(X, Y, Z)

# ╔═╡ 6cac6c4b-a0d7-4cae-bc2d-1cf9e0d775a2
md"#### Example 1: Build a Heisenberg Hamiltonian"

# ╔═╡ 1b35ea7b-c831-4fce-8109-3cd35a80214a
heisenberg(n) = sum([sum([put(n, i=>G) * put(n, i+1=>G) for G in [X, Y, Z]]) for i=1:n-1])

# ╔═╡ a4b4761b-0686-42bd-9f04-d6083b6064e4
heisenberg(3)

# ╔═╡ f639e1b4-e257-4218-8727-292eadedd7bf
mat(heisenberg(20))

# ╔═╡ 6aec6b4e-f904-4e3d-9279-1507e0d39eec
md"## Measurement"

# ╔═╡ e46611a9-2552-463b-b6f2-b3ad8abd28fc
md"readout results -> measure"

# ╔═╡ 1b9ff9b0-fd59-41fc-98b0-650eddf52ab2
md"""
```math
\langle O\rangle = \langle\psi|O|\psi\rangle
```
"""

# ╔═╡ 6294a05c-46b1-47e4-9f3c-e741ab191a09
obs = kron(Z, Z)

# ╔═╡ 22175f6d-b984-4bfc-b776-5f6da6742e75
Yao.ishermitian(obs)

# ╔═╡ 0139dbd2-486e-4b58-b656-6b5e06864cd1
expectation_value = expect(obs, ψt)  # not physically implementable

# ╔═╡ f02617e1-8856-4369-a33d-c73e555435b2
md"Measure is a bit complicated"

# ╔═╡ 8a94551d-c82a-4db0-aaf0-3c5285ab0695
let
	ψt2 = copy(ψt)
	measure!(ψt2), statevec(ψt2)
end

# ╔═╡ 9eb46218-30c9-4109-9b32-9382fdee1081
let
	ψt2 = copy(ψt)
	measure!(obs, ψt2), statevec(ψt2)
end

# ╔═╡ 1ef52a16-16f9-40c6-b61c-3a00dfa3ce9f
md"# Quantum circuit simulation"

# ╔═╡ f1bc3f79-5c82-4e4f-b970-81ea33c4493b
md"""
1. quantum state -> quantum bits
2. time evolution ->  applying the circuit
3. measure  ->  readout
"""

# ╔═╡ a74105a9-0a55-419b-895b-f3c2c831bffb
a_random_circuit = chain(dispatch!(variational_circuit(5, 2), :random), Yao.Measure());

# ╔═╡ 1465cd17-94f7-44b9-9179-4d7b39ae8c3a
vizcircuit(a_random_circuit; starting_texts=["0" for i=1:5],  ending_texts=["?" for i=1:5])

# ╔═╡ c09975ac-4387-451d-98e7-800b8b14760c
statevec(zero_state(5))

# ╔═╡ 2bf0685a-42ff-4d73-a009-6b1523ce7f4b
function plot_amplitude(reg::ArrayReg)
	Plots.bar(0:1<<nqubits(reg)-1, reg |> probs, label="probability")
end;

# ╔═╡ 27ee8d46-710c-48fe-981c-6ab282843741
function plot_measured(measure_results::Vector{<:BitStr{N}}) where N
	bar(0:1<<N-1, [count(b->Int(b)==i, measure_results) for i=0:1<<N-1], label="counting")
end;

# ╔═╡ 1224f1bc-9835-42eb-9b9c-4e5628575a4c
reg0 = zero_state(5)

# ╔═╡ 38cf5a1a-e995-4949-bd4e-deba6500cf00
plot_amplitude(reg0)

# ╔═╡ c2e846cb-9f80-4ea1-8295-eb9dadd7b878
md"Measuring on computational basis ``\{P_1^{(n)}, P_1^{(n-1)}, \ldots P_1^{(1)}\}``"

# ╔═╡ 480e476f-c62e-49e6-b022-4f10c5b42e57
measure!(reg0)

# ╔═╡ 922f19c7-4d2f-48a7-9fff-4c63acb868ba
let
	reg0 = zero_state(5)
	plot_amplitude(reg0 |> a_random_circuit)
	Plots.annotate!(30, 0.5, "$(bitstring(Int(a_random_circuit[end].results))[end-4:end])")
end

# ╔═╡ 812aa02f-e8b0-4226-b2ef-5260884b2f09
function run_and_get_measure_results(circuit)
	m = circuit[end]
 	zero_state(nqubits(circuit)) |> circuit
	m.results
end

# ╔═╡ fe69d633-1da8-4db5-9e30-eb19c4ec392f
 measure_results = [run_and_get_measure_results(a_random_circuit) for i=1:1000]

# ╔═╡ 23a5b68f-bc93-43a6-9072-7d929f5f62dc
plot_measured(measure_results)

# ╔═╡ 847051a1-48e9-4852-8f6d-9a5f604e9142
md"## Deutsch–Jozsa algorithm"

# ╔═╡ 85228d6e-731d-465d-a000-672e7fff8aff
md"![](https://upload.wikimedia.org/wikipedia/commons/thumb/b/b5/Deutsch-Jozsa-algorithm-quantum-circuit.png/600px-Deutsch-Jozsa-algorithm-quantum-circuit.png)"

# ╔═╡ f66bd812-cfd0-4ff2-be21-71c38cb8ab79
md"Case 1, ``f(x)`` is a constant.
"

# ╔═╡ 6d6b55d4-4ddb-4ef9-b2b1-154021c01273
md"Case 2, ``f(x)`` is uniform."

# ╔═╡ 5467dafc-265a-4848-a898-364124a6e9b4
mat(chain(H, H))

# ╔═╡ 1e01ef7f-63d0-4b0d-8455-7b524206f15f
statevec(zero_state(1) |> X |> H)

# ╔═╡ 1ca545b9-6f88-4251-aa8a-d5ad0236c77e
statevec(zero_state(1) |> X |> H |> X)

# ╔═╡ 528bfa77-d6da-4700-a72c-5afffdbe6139
zero_state(4) |> repeat(4, H, 1:4) |> statevec

# ╔═╡ 1571a8ae-aae5-4af8-9b9c-6fe9e832a39e
md"""
```math
\begin{align}
|\psi\rangle &= \sum_x|x\rangle \otimes \left(X^{f(x)}\frac{|0\rangle + |1\rangle}{\sqrt{2}}\right)\\
 &= \sum_{x\in\{x|f(x)=0\}}|x\rangle \otimes \frac{|0\rangle + |1\rangle}{\sqrt{2}} - \sum_{x\in\{x|f(x)=1\}}|x\rangle \otimes \frac{|0\rangle + |1\rangle}{\sqrt{2}}\\
&= \underbrace{\sum_{x\in\{x|f(x)=0\}}|x\rangle - \sum_{x\in\{x|f(x)=1\}}|x\rangle}_{\text{has zero overlap with $(H |0\rangle)^{\otimes n}$}} \otimes \frac{|0\rangle + |1\rangle}{\sqrt{2}}
\end{align}
```
"""

# ╔═╡ c60ea349-6b96-4f70-8369-a92379338bbd
deutsch_jozsa(circuit::AbstractBlock{N}, inputs, output) where N = chain(
	put(N+1, N+1=>X), 
	repeat(N+1, H, [inputs..., N+1]),
	subroutine(N+1, circuit, 1:N),
	control(N+1, output, N+1=>X), 
	subroutine(N+1, circuit', 1:N),
	repeat(N+1, H, inputs),
	Yao.Measure(N+1, locs=inputs)
)

# ╔═╡ d8e9ae38-1569-497a-8506-a3e059ffc7ab
let
	c = deutsch_jozsa(constantf(4), 1:4, 5)
	vizcircuit(c; scale=0.5)
end

# ╔═╡ 195ccd7e-f1ea-4704-9392-ffe596a756f9
let
	addc, outputs = add_circuit(2)
	c = deutsch_jozsa(addc, 1:4, outputs[1])
	zero_state(nqubits(c)) |> c
	#vizcircuit(c; scale=0.3)
	c[end].results
end

# ╔═╡ f1814dfa-588d-49b8-a68e-bda185f3712f
let
	c = deutsch_jozsa(constantf(4), 1:4, 5)
	zero_state(nqubits(c)) |> c
	#vizcircuit(c)
	c[end].results
end

# ╔═╡ 640a099a-2af9-4493-a103-9a371c1e32ba
let
	c = deutsch_jozsa(constantf(4), 1:4, 5)[1:end-1]  # we do not want the measure
	expect(repeat(nqubits(c), P0, 1:4), zero_state(nqubits(c)) |> c)
end

# ╔═╡ 56aa1520-8289-4082-adc7-c12c2b86875b
let
	fc, outputs = add_circuit(2)
	n = nqubits(fc)
	c = deutsch_jozsa(fc, 1:4, outputs[1])[1:end-1]  # we do not want the measure
	#vizcircuit(c; scale=0.3)
	expect(repeat(n+1, P0, 1:4), zero_state(n+1) |> c)
end

# ╔═╡ 01f0cd17-ed40-4b7a-93a6-28f80c1367ea
md"## Quantum is not only probability"

# ╔═╡ 6271ddf9-dea1-4b69-baed-ce4387fa91d7
zero_state(4) |> repeat(4, H, 1:4) |> probs

# ╔═╡ 84d4d1fd-0921-48ea-b292-fa2ccce58e7d
zero_state(4) |> repeat(4, H, 1:4) |> put(4, 1=>X) |> probs

# ╔═╡ 952bd5bf-c7d9-483c-b8e7-a6e371808548
let
	reg1 = zero_state(4) |> repeat(4, H, 1:4) |> put(4, 1=>Z)
	reg2 = zero_state(4) |> repeat(4, H, 1:4)
	reg1' * reg2
end

# ╔═╡ c0e0ee8e-8580-424a-8d9a-33f2c8ac47c5
md"""
![](https://cdn1.byjus.com/wp-content/uploads/2020/11/Young_s-Double-Slit-Experiment-9.png)
"""

# ╔═╡ 840e56f4-2fcf-4136-bdb7-7feba88f25da
md"# Tensor network based quantum simulation"

# ╔═╡ 09137c47-1bd0-481b-8b60-ca2abff64afc
md"""## Exponentially large vector is not the only solution to quantum simulation
```math
|\psi_1\rangle = v_1 \otimes v_2  \otimes \ldots \otimes v_n
```
```math
v_i = \left(\begin{matrix}a_i\\ b_i\end{matrix}\right)
```
"""

# ╔═╡ d3870b24-b9f7-44e9-8d9f-0408978df4e7
md"Why not evaluate it directly?"

# ╔═╡ c2811274-f59c-4021-b939-24ab695c24ba
md"""
```math
|\psi_2\rangle = w_1 \otimes w_2  \otimes \ldots \otimes w_n
```
"""

# ╔═╡ 13d5d5ef-28c2-40e3-b686-cf2db0ace884
md"""
```math
\langle \psi_2|\psi_1\rangle = \sum_{s \in \{0, 1\}^{n}} \prod_{i=1}^n (v_i)_{s_i} (w_i)_{s_i}
```
"""

# ╔═╡ a55ded4b-1ecc-480f-97f9-bec5d7b475f6
md"There are exponential number of terms!"

# ╔═╡ 463cf35c-d9c7-437b-858e-8fde0d37ae58
md"""
```math
\langle \psi_2|\psi_1\rangle = \prod_{i=1}^n \sum_{s_i \in \{0, 1\}} (v_i)_{s_i} (w_i)_{s_i}
```
"""

# ╔═╡ 2fad56df-cb36-4eb4-a853-ac3b22de86b0
md"Only linear number of terms"

# ╔═╡ 0854f2c0-ecb2-4e21-bebc-f6e480ca3aa6
canvas() do
	Compose.set_default_graphic_size(5cm, 5cm)
	nb = nodestyle(:default)
	line = bondstyle(:default)
	tb = textstyle(:default)
	for i=1:5
		n1 = nb >> (0.2, i*0.2-0.1)
		tb >> ((0.1, i*0.2-0.1), "v$i")
		
		n2 = nb >> (0.8, i*0.2-0.1)
		tb >> ((0.9, i*0.2-0.1), "w$i")
		line >> (n1, n2)
		tb >> ((0.5, i*0.2-0.2+0.05), "s$i")
	end
end

# ╔═╡ 5c0d6fa7-bfc7-4d34-9b2c-93b49df28b17
md"## A quantum circuit is a sum product network!"

# ╔═╡ 812b1a90-83fb-4c4f-8ac5-e8d00b13972c
md"Sum product notation for quantum computing 
```math
\langle \psi_2|U|\psi_{1}\rangle
```
where ``U`` is your circuit.
"

# ╔═╡ 363e5f54-1659-4953-8bc5-3d05911126f3
canvas() do
	Compose.set_default_graphic_size(5cm, 5cm)
	nb = nodestyle(:default)
	line = bondstyle(:default)
	tb = textstyle(:default)
	for i=1:5
		n1 = nb >> (0.1, i*0.2-0.1)
		n2 = nb >> (0.9, i*0.2-0.1)
		line >> (n1, n2)
	end
	YaoPlots.CircuitStyles.G() >> (0.5, 0.4)
	YaoPlots.CircuitStyles.TEXT() >> ((0.5, 0.4), "U")
	for loct in [((0.5, 0.05), "a"), ((0.2, 0.25), "b"), ((0.8, 0.25), "β"),
			((0.2, 0.45), "c"), ((0.8, 0.45), "γ"), ((0.5, 0.65), "d"), ((0.5, 0.85), "e")]
		tb >> loct
	end
end

# ╔═╡ 5d8943a0-554e-47ba-ba9c-1ccf0387d844
code = ein"a,b,c,d,e,βγbc,a,β,γ,d,e->"

# ╔═╡ d107918e-7ee1-4f1c-b11c-7d3a8821f5ae
vs, ws = [rand_state(1) for i=1:5], [rand_state(1) for i=1:5]

# ╔═╡ bc40176b-08ee-478c-96c6-1b15b4630c78
U = control(2, 1, 2=>X)

# ╔═╡ c99cc454-8475-49fa-9a27-e456b93926a6
ein"a,b,c,d,e,βγbc,a,β,γ,d,e->"(statevec.(vs)..., reshape(Matrix(U), 2, 2, 2, 2), conj.(statevec.(ws))...)

# ╔═╡ 95c3af25-5771-453b-aede-1878932d0e6e
let
	reg = join(vs[end:-1:1]...)
	reg2 = join(ws[end:-1:1]...)
	reg2' * (reg |> put(5, (2,3)=>U))
end

# ╔═╡ c59f5370-98c9-49fc-abd0-b9a0f61a9fe2
timespace_complexity(code, uniformsize(code, 2))

# ╔═╡ 9cc7d3ef-a26f-469d-8295-b248936a31e2
timespace_complexity(optimize_code(code, uniformsize(code,2), TreeSA()), uniformsize(code, 2))

# ╔═╡ 62393bfb-5597-40cd-a5e4-c4b10a61d1cf
md"""
#### The evaluation order of sums and products matters!
* Simulating quantum computation by contracting tensor networks, Igor L. Markov, Yaoyun Shi (arXiv:quant-ph/0511069)
* Solving the sampling problem of the Sycamore quantum supremacy circuits, Feng Pan, Keyang Chen, Pan Zhang (arXiv:2111.03011)
* Limitations of Linear Cross-Entropy as a Measure for Quantum Advantage, Xun Gao, Marcin Kalinowski, Chi-Ning Chou, Mikhail D. Lukin, Boaz Barak, Soonwon Choi (arXiv:2112.01657)
"""

# ╔═╡ ed4c058f-029b-4162-8d8b-5c323418e8ad
md"## When tensor network is good."

# ╔═╡ a06024c9-49f5-44ab-a3e5-9c06ee68f34c
md"![](https://user-images.githubusercontent.com/6257240/146654677-09b902f3-a32a-4695-97f2-598199419142.png)

Source: Variational Quantum Eigensolver with Fewer Qubits
Jin-Guo Liu, Yi-Hong Zhang, Yuan Wan, Lei Wang (arXiv:2112.01657)
"

# ╔═╡ 7a39257d-b864-460d-bccd-eed9933c4b7c
function mps_like_circuit(n)
	chain(repeat(n, H, 1:n), [control(n, i, i+1=>Z) for i=1:n-1]...)
end

# ╔═╡ 19c37ef5-fe1d-4061-a5ad-f4b6bd9ecba4
mps_like_circuit(4) |> vizcircuit

# ╔═╡ d04690ca-390b-462b-8257-e9ebe01b3fd0
md"""
# Summary

* The magic of quantum computing is related to the $(highlight("Interference")) effect,
* Quantum simulation requires $(highlight("exponential")) resources,
-------------------
* A quantum circuit is a $(highlight("sum-product")) network,
* The contraction order of a sum-product network matters,
* sum-product network is good at handling $(highlight("low entangled")) states.
-------------------
* Two good packages: Yao.jl and OMEinsum.jl (star us on Github!).
### How to find me
* Github issue (Yao.jl, OMEinsum.jl)
* Julia slack (channel: #yao-dev, or @GiggleLiu)
"""

# ╔═╡ e8201e38-73dc-45d4-b93d-20f063d2ac3a
md"## Quantum Compiling"

# ╔═╡ aa45acea-e6f7-49fb-a5b6-3ea0f2c1530c
begin
	function decompose_toffoli(x::ControlBlock{N,XGate,2,1}) where N
		if x.ctrl_config == (1, 1)
			i, j, k = x.ctrl_locs[1], x.ctrl_locs[2], x.locs[1]
			return chain(put(N, k=>H), control(N,j,k=>X), put(N, k=>ConstGate.T'), control(N, i, k=>X), put(N, k=>ConstGate.T), control(N, j, k=>X), put(N, k=>ConstGate.T'), control(N, i, k=>X), put(N, k=>ConstGate.T), put(N, j=>ConstGate.T), control(N, i, j=>X), put(N, k=>H), put(N, i=>ConstGate.T), put(N, j=>ConstGate.T'), control(N, i, j=>X))
		else
			return x
		end
	end
	decompose_toffoli(x::AbstractBlock) = chsubblocks(x, decompose_toffoli.(subblocks(x)))
end

# ╔═╡ ac82ea9c-be86-4501-95db-04f7eef26bea
cluster_code = let
	c = decompose_toffoli(mps_like_circuit(10))
	code, xs = yao2einsum(c; initial_state=Dict([i=>zero_state(1) for i=1:nqubits(c)]), final_state=Dict([i=>zero_state(1) for i=1:nqubits(c)]))
	optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
end

# ╔═╡ b931b25c-508e-4e64-bc3a-ff993643e165
timespace_complexity(cluster_code, uniformsize(cluster_code, 2))

# ╔═╡ fdf13c95-6d95-4220-9759-5fea44eea274
optcode = let
	c = decompose_toffoli(add_circuit(4)[1])
	code, xs = yao2einsum(c; initial_state=Dict([i=>rand(0:1) for i=1:nqubits(c)]), final_state=Dict([i=>rand(0:1) for i=1:nqubits(c)]))
	optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
end

# ╔═╡ 669bed15-1c57-401b-bcfb-48b5a1520a58
timespace_complexity(optcode, uniformsize(optcode, 2))

# ╔═╡ 1e2848d6-252b-4db5-9864-1f1886da9998
vizcircuit(decompose_toffoli(control(3, (1,2), 3=>X)); w_depth=0.7, scale=0.6)

# ╔═╡ e2c3c7a8-0349-4f4c-87a7-3fbe8fc2fe46
vizcircuit(decompose_toffoli(add_circuit(4)[1]); scale=0.3)

# ╔═╡ 0e11d98d-79c6-444d-8df1-357e8218f233
gatecount(decompose_toffoli(add_circuit(4)[1]))

# ╔═╡ 800262e6-5c78-4229-82d5-40de619d3013
md"#### Example 2: Build a quantum fourier transformation circuit"

# ╔═╡ a6db5c0e-ea51-4071-b3fc-146ef90270aa
cphase(i, j) = control(i, j=> shift(2π/(2^(i-j+1))));

# ╔═╡ ecf006d7-335c-48de-b629-b561dc08b334
hcphases(n, i) = chain(n, i==j ? put(i=>H) : cphase(j, i) for j in i:n);

# ╔═╡ 22fba19f-d337-48cc-a311-e44d03a0c050
qft_circuit(n::Int) = chain(n, hcphases(n, i) for i = 1:n)

# ╔═╡ c41dc6d8-07ff-42fb-a54a-d1a1cf1cc223
vizcircuit(qft_circuit(5))

# ╔═╡ 70e85333-26b7-40b8-8fe1-7fb470b5f7b3
let
	vizcircuit(qft_circuit(5))
end

# ╔═╡ 28f657cd-8ae7-4572-94fa-813797b59b25
import Cairo

# ╔═╡ Cell order:
# ╠═4a96f5c9-37b4-4a8a-a6bd-8a4b4440eb49
# ╠═2a145cba-26b0-43bd-9ab2-13818d246eae
# ╠═e11c26c0-e534-45fc-bb1c-c0f2ce4016db
# ╠═a8b06352-5265-4b89-bd82-b31f3cdac391
# ╟─e6306a69-bd6a-4c01-9c6e-1cb523668019
# ╟─86b15fb0-f112-4689-8663-8cc6c0a8fb2a
# ╟─0f099c85-f039-477e-a70d-a3801cbb2656
# ╟─48854a73-4896-4542-9ad4-15ae87418f1d
# ╟─c1d40103-1710-4221-b414-0958c13fb95f
# ╟─4c8a4cea-63f4-49f0-82b9-dfa608be46bf
# ╟─eb75810a-a746-4e4a-889f-c87c8d1e153f
# ╟─330356f0-d7d7-425f-b16c-d1e9f97494b5
# ╟─b83f7675-68fa-44b5-8681-c85984eeb877
# ╟─b7d90c25-2f66-4ace-8d52-841ea376b3f9
# ╟─2e9c5b7f-14dc-4ea1-a286-6a68c3e96b3e
# ╟─59c24944-9a8d-4a7f-8b0b-c08a000cd655
# ╟─2c6d024a-3187-4976-af27-393af8826a2d
# ╟─165b369a-3d0a-4ae7-99f1-de0297f93707
# ╟─c841f2e2-b907-4f74-be57-968ca339bec4
# ╟─e1b6b7a9-2d26-4f43-a0bd-54f7ad22a5b3
# ╟─16febe07-7f52-4b74-b9d6-e20fd8b05ab3
# ╟─d631887a-a222-402d-94bb-ecc6db6bea56
# ╟─42f3e654-1836-491a-8119-b03b93822f45
# ╟─cc0d4b1c-07b6-42ed-b1e2-80581a02ee6f
# ╟─c5694b66-d023-42c1-ae62-b7218ba8ebe7
# ╟─c5812ca6-c4ca-4211-9d2d-df498fd7a2da
# ╟─6b4e2611-5a5a-4646-aa8b-fec0d84d240c
# ╟─ad47b86e-d8da-401d-a1a2-1139664adee5
# ╟─ed0dc5f3-0da3-4d4f-9e3b-d4f471aa2e03
# ╟─ec1613ad-482e-4c68-b130-b51725f8d94e
# ╟─bde5a21b-7db3-40c2-8c00-54182987dcfd
# ╟─07f36b6e-a98d-4bee-9091-b1a717dab0e8
# ╟─a4bd7a89-7dda-4652-a30d-3e2b02079b1e
# ╟─91f767ea-5afa-4e04-9101-d1c9b45a5b1e
# ╟─e188093d-1ef1-4e70-966b-77cc0761a801
# ╟─f6533ed1-b419-43a7-b03e-3a538d2e46da
# ╟─cac05229-df74-4679-89ab-e65003b7d773
# ╟─4518fcfa-9f04-4c2c-82b8-719033f20ac9
# ╟─bded6dc0-bce0-404e-a1f2-36712f26f2d6
# ╟─beed43ce-b24c-4edd-9da2-30bdd72c9411
# ╟─0c267c28-967b-4268-8372-00d7b48b7b8e
# ╟─7440555d-685f-4447-a44a-0463b37bd43c
# ╟─fed37413-c341-4dcf-a637-404d2f186b2b
# ╟─69d1fc34-3e02-4c4b-83a7-6d05a0d69405
# ╟─67ad78c7-930a-48b1-b8d3-d9e549ee7379
# ╟─523fdddf-e76c-46f0-bda0-4d3e19c7a82d
# ╟─f311befa-ba09-42e6-ac9b-b59450162ebd
# ╟─816ebc7b-43ec-47e2-a24a-86ac70dd6afe
# ╟─893baeda-d5b0-43c1-9018-8d18c06486ca
# ╟─9a5227f6-0db8-43c1-9b79-80897faa86d0
# ╟─cee8db91-7e84-4728-aef6-ea861c62ff96
# ╟─38194214-3bf1-4229-9fe8-37282b30a5ad
# ╟─4591b3a0-ff3f-46c8-8c8d-ef7c6ba76cc5
# ╟─7c1cd22b-c478-4dcf-9828-ae6aad3df95a
# ╟─ae12828c-350a-4251-96cf-5ac2dc0fc0ee
# ╟─6cc597c9-e370-48ea-a5c5-422904d4d8a0
# ╟─b2bf15bb-1350-4304-a932-c87b09558115
# ╟─0e63c1eb-fea6-411a-bae6-0bc90dee6bc7
# ╟─0b42d7bb-b9cf-42ef-aad7-5fc5d5918be3
# ╟─622fdaac-87dd-4c22-85b9-470510566480
# ╟─26e5d619-7119-49bc-8907-17ae0db424f5
# ╟─a7d85b82-a705-4f3b-a371-06a87071335d
# ╟─86b398d0-844d-4d06-9ede-79d618502cce
# ╟─e2e4ec60-2737-4560-89b1-1e14a35044e8
# ╟─aea490af-0bdd-4930-9ad2-7d9a13e08c46
# ╟─7f84eed8-edd8-4f8e-a2b1-3ad862285934
# ╟─c0078012-3d81-4584-a050-9a58802d08a9
# ╟─ea458fa2-1f9f-46e1-88da-942034d0fa73
# ╟─3381001d-1120-4b88-ac01-5ca861f0a9be
# ╟─7524029c-2b7e-465e-9827-c993d6cdd34a
# ╟─22cd8a36-c76f-4a8c-a8a2-e1924136012a
# ╟─051548ff-dd39-4d55-ac53-e8e2bacec68e
# ╟─e7b2e327-952e-4062-9193-19653eaee19c
# ╟─cb2c496c-e9dc-4666-b41b-d97cca377047
# ╟─78cb38db-b0b6-49f7-bcaf-1a2df58b8150
# ╟─334db854-7805-4548-a414-ad7d215fe387
# ╟─40956c0e-3833-4397-a4a8-c46e1890ef39
# ╟─8bebf527-d197-48b0-9778-f71b0a5dad77
# ╟─ca094f81-bef3-40eb-9279-086fe7eb506a
# ╟─0840e7ea-aa63-44df-a788-ad18ac842006
# ╟─e5004508-2b9a-4acd-91ba-ad21fc7d1b16
# ╟─5c03fc94-ce6d-490b-a3ee-911692877e6a
# ╟─2da83f06-e639-4f75-a162-cb561a8207ca
# ╟─960b0cdc-a769-4f03-8842-4e9d86e67cef
# ╟─408f534d-5553-4a95-8d14-cabf36f8aff5
# ╟─41244f02-796b-4afb-9825-0b21d216b78b
# ╟─b8069c02-d422-4162-ad9c-c4b94870db6d
# ╟─57a3616e-49af-40b7-a000-a4ecc81af84e
# ╟─0f655562-0c18-43a8-83f0-4d11b7eeb068
# ╟─774def48-8f9f-425a-aa47-ffe396692247
# ╟─8017edf3-05b2-4cee-8a4e-90331a681037
# ╟─d11e3541-6e03-41d7-a04b-fcce244058f7
# ╟─33bd9c30-8570-4580-94a7-adbfec13bedc
# ╟─40a3bc85-be6d-4d1b-8376-da5a33c7057b
# ╟─52d14933-99d4-4ee6-9eff-7be4a5722334
# ╟─e5d6f98c-aa84-48cf-bbca-cd72333d67f7
# ╠═40563837-b1e9-4df9-9cc7-8841d0068973
# ╟─99b9d334-40ed-4636-b299-1afe40402a12
# ╟─400f697d-fbdb-4067-881f-e004f50bbc0f
# ╠═b3daf0a3-4d1e-4474-bb99-948c718059ce
# ╠═30c16e7d-417c-4098-a879-0914d711ff56
# ╟─6a719cf5-d8a4-4369-95e7-c1dfbc5cebf8
# ╠═6f64057b-130d-44c1-8957-e4c7c56f7cb4
# ╟─417b6e83-4320-40eb-a91c-9e6cf5f36be8
# ╠═cd39c31a-6fc9-489f-8045-cd3eca229a76
# ╠═8d1c023b-4ee6-4201-99f6-0245cee70cda
# ╠═7278c646-3c7b-4b18-a989-cb98abefd5cc
# ╟─30100576-9fa2-4825-bfa6-35ae6cb920bb
# ╟─217e42d4-d226-42b2-bbd7-de85bc467ec4
# ╟─0e1e5385-db28-474f-a14a-cf3fd1efcd94
# ╟─681cdc3b-348f-45c4-9031-6298e81d995d
# ╠═b3813f2d-4e9d-4789-a653-d28fd1b98eda
# ╠═a9d19c74-1b83-4211-9c83-95ad1059d432
# ╠═48427378-face-48eb-b635-f3f460a51e26
# ╟─45d19896-ae93-4f9a-b1fe-bb143c7ff580
# ╟─bf983162-66e3-40e0-afce-7f5b7728c490
# ╠═bf357f44-54e2-49c0-b0b6-263abb1dbcb0
# ╠═913bc173-8cfd-4680-ad13-36ea806acecb
# ╠═675b3398-01d7-4949-bb0e-7cdf9b805c69
# ╠═09df3c8e-3d25-4ba2-bf19-9947224dfca4
# ╠═45ace908-90df-488c-9258-459c0344baae
# ╠═de21a067-4587-461d-94bb-a295b6768f96
# ╠═d85a2028-9b6f-4d6c-ba61-fecef3fd4876
# ╠═4a7e56ac-a69f-4b1f-b585-3a92f43e42e8
# ╟─f4f18b91-86f6-460c-b5de-4c52931a4098
# ╟─bfb9fcb5-2555-490a-9dcc-48cc91f8ab3c
# ╠═2e17662d-24ed-40c9-93b0-4cef526c3a75
# ╠═21b11fac-5efb-4fd5-a1b0-e684d215a46c
# ╟─3fa73e8a-fe46-475f-8346-e54f52c144f4
# ╠═31c1c4fc-7476-4066-8462-bf18f8d69966
# ╟─81a55867-4158-4371-a12c-b08e3d64c0cb
# ╠═7f15cbb0-09be-46ee-852c-e43dde9bc4f5
# ╟─cc4e7d7b-0596-4bc0-b23e-f6861fcd5260
# ╠═8b03f2f1-f38c-4d5e-b778-0c3b4aaf910d
# ╠═a2129e65-7e73-4b42-9924-e88d60893ed2
# ╟─4dd7f18c-2c0f-4405-bbfe-706699ea8958
# ╟─faf2b47d-ab3b-4087-8cda-a564303c3bc9
# ╠═8b1606b5-78d8-4020-8597-1cf9f2364a2b
# ╠═04693ee8-d0b2-4ffa-80b2-40b9149eceea
# ╠═8be8d92a-d882-4e77-a3c3-9846babeda2b
# ╟─597bff8b-f7fa-4fc8-bc59-f06c934f211e
# ╟─421abcef-a1be-4297-a159-adcd92332d6b
# ╠═3388350c-fd8b-4cae-bf34-71446bd747ca
# ╠═99847765-fd49-4c9d-914b-5d47eac89490
# ╠═fbd1ea6e-16da-4289-b0e3-d8e8a5a55541
# ╠═5dafbbbd-d6f7-423c-8a46-41b1596f83a3
# ╠═40889c09-1a10-4829-9a4e-5dc46d9045ef
# ╠═984db7f2-30e0-421d-9f27-874553b6fae2
# ╟─44d6b7ba-43e6-43f7-bbed-a727edd22e0c
# ╟─f35b9789-788b-450d-af8e-472412efce77
# ╠═34a0b380-ac20-4d88-b498-969745a0df7c
# ╠═3aaa7e49-9507-4384-ac7e-f4d486127811
# ╟─ac2f9c18-f4bc-4077-9ad1-1f4f0d75a7b1
# ╠═47da74c7-0f24-4dd5-b6f7-88c22a81de5a
# ╠═e06b23aa-7fae-4495-9b69-43cea3682a28
# ╠═915185d6-3b36-437e-bac0-22cd5503bba2
# ╟─67d9d529-aeec-4ede-a0f7-85dc7f4fba5a
# ╟─0a6d5a41-f42b-42e7-adf7-93a40f67c9f3
# ╠═dc46cbbe-0eb5-41f9-ad5f-f8a7fba5ff9c
# ╠═0a7cab72-a50b-4d58-877a-4a9d939b198f
# ╠═6b2d8bb9-5a4f-4699-ac92-bc8a4003918f
# ╟─17a574c9-664a-4cba-8571-60a2551f37f7
# ╟─3a89d593-ae7e-454c-bb18-0baa4cc06a25
# ╠═4559f70b-5a9b-478d-ba9a-3a35f6fc9716
# ╟─6cac6c4b-a0d7-4cae-bc2d-1cf9e0d775a2
# ╠═1b35ea7b-c831-4fce-8109-3cd35a80214a
# ╠═a4b4761b-0686-42bd-9f04-d6083b6064e4
# ╠═f639e1b4-e257-4218-8727-292eadedd7bf
# ╟─6aec6b4e-f904-4e3d-9279-1507e0d39eec
# ╟─e46611a9-2552-463b-b6f2-b3ad8abd28fc
# ╟─1b9ff9b0-fd59-41fc-98b0-650eddf52ab2
# ╠═6294a05c-46b1-47e4-9f3c-e741ab191a09
# ╠═22175f6d-b984-4bfc-b776-5f6da6742e75
# ╠═0139dbd2-486e-4b58-b656-6b5e06864cd1
# ╟─f02617e1-8856-4369-a33d-c73e555435b2
# ╠═8a94551d-c82a-4db0-aaf0-3c5285ab0695
# ╠═9eb46218-30c9-4109-9b32-9382fdee1081
# ╟─1ef52a16-16f9-40c6-b61c-3a00dfa3ce9f
# ╟─f1bc3f79-5c82-4e4f-b970-81ea33c4493b
# ╠═0f8d63ac-f677-4889-b033-2a93f62be700
# ╠═a74105a9-0a55-419b-895b-f3c2c831bffb
# ╠═1465cd17-94f7-44b9-9179-4d7b39ae8c3a
# ╠═c09975ac-4387-451d-98e7-800b8b14760c
# ╟─2bf0685a-42ff-4d73-a009-6b1523ce7f4b
# ╟─27ee8d46-710c-48fe-981c-6ab282843741
# ╠═1224f1bc-9835-42eb-9b9c-4e5628575a4c
# ╠═38cf5a1a-e995-4949-bd4e-deba6500cf00
# ╟─c2e846cb-9f80-4ea1-8295-eb9dadd7b878
# ╠═480e476f-c62e-49e6-b022-4f10c5b42e57
# ╟─922f19c7-4d2f-48a7-9fff-4c63acb868ba
# ╠═812aa02f-e8b0-4226-b2ef-5260884b2f09
# ╠═fe69d633-1da8-4db5-9e30-eb19c4ec392f
# ╠═23a5b68f-bc93-43a6-9072-7d929f5f62dc
# ╠═342685aa-5159-11ec-13fd-fb8954106bca
# ╟─847051a1-48e9-4852-8f6d-9a5f604e9142
# ╟─85228d6e-731d-465d-a000-672e7fff8aff
# ╟─f66bd812-cfd0-4ff2-be21-71c38cb8ab79
# ╟─6d6b55d4-4ddb-4ef9-b2b1-154021c01273
# ╠═5467dafc-265a-4848-a898-364124a6e9b4
# ╠═1e01ef7f-63d0-4b0d-8455-7b524206f15f
# ╠═1ca545b9-6f88-4251-aa8a-d5ad0236c77e
# ╠═528bfa77-d6da-4700-a72c-5afffdbe6139
# ╟─1571a8ae-aae5-4af8-9b9c-6fe9e832a39e
# ╠═c60ea349-6b96-4f70-8369-a92379338bbd
# ╠═d8e9ae38-1569-497a-8506-a3e059ffc7ab
# ╠═195ccd7e-f1ea-4704-9392-ffe596a756f9
# ╠═f1814dfa-588d-49b8-a68e-bda185f3712f
# ╠═5fbddde1-3594-43c4-9f56-c0ae258d926f
# ╠═640a099a-2af9-4493-a103-9a371c1e32ba
# ╠═56aa1520-8289-4082-adc7-c12c2b86875b
# ╟─01f0cd17-ed40-4b7a-93a6-28f80c1367ea
# ╠═6271ddf9-dea1-4b69-baed-ce4387fa91d7
# ╠═84d4d1fd-0921-48ea-b292-fa2ccce58e7d
# ╠═952bd5bf-c7d9-483c-b8e7-a6e371808548
# ╟─c0e0ee8e-8580-424a-8d9a-33f2c8ac47c5
# ╟─840e56f4-2fcf-4136-bdb7-7feba88f25da
# ╟─09137c47-1bd0-481b-8b60-ca2abff64afc
# ╟─d3870b24-b9f7-44e9-8d9f-0408978df4e7
# ╟─c2811274-f59c-4021-b939-24ab695c24ba
# ╟─13d5d5ef-28c2-40e3-b686-cf2db0ace884
# ╟─a55ded4b-1ecc-480f-97f9-bec5d7b475f6
# ╟─463cf35c-d9c7-437b-858e-8fde0d37ae58
# ╟─2fad56df-cb36-4eb4-a853-ac3b22de86b0
# ╟─0854f2c0-ecb2-4e21-bebc-f6e480ca3aa6
# ╟─5c0d6fa7-bfc7-4d34-9b2c-93b49df28b17
# ╟─812b1a90-83fb-4c4f-8ac5-e8d00b13972c
# ╟─363e5f54-1659-4953-8bc5-3d05911126f3
# ╠═5d8943a0-554e-47ba-ba9c-1ccf0387d844
# ╠═d107918e-7ee1-4f1c-b11c-7d3a8821f5ae
# ╠═bc40176b-08ee-478c-96c6-1b15b4630c78
# ╠═c99cc454-8475-49fa-9a27-e456b93926a6
# ╠═95c3af25-5771-453b-aede-1878932d0e6e
# ╠═c59f5370-98c9-49fc-abd0-b9a0f61a9fe2
# ╠═9cc7d3ef-a26f-469d-8295-b248936a31e2
# ╟─62393bfb-5597-40cd-a5e4-c4b10a61d1cf
# ╟─ed4c058f-029b-4162-8d8b-5c323418e8ad
# ╟─a06024c9-49f5-44ab-a3e5-9c06ee68f34c
# ╠═6f131f36-0b4e-4570-8527-620297fae48e
# ╠═7a39257d-b864-460d-bccd-eed9933c4b7c
# ╠═19c37ef5-fe1d-4061-a5ad-f4b6bd9ecba4
# ╠═ac82ea9c-be86-4501-95db-04f7eef26bea
# ╠═b931b25c-508e-4e64-bc3a-ff993643e165
# ╠═fdf13c95-6d95-4220-9759-5fea44eea274
# ╠═669bed15-1c57-401b-bcfb-48b5a1520a58
# ╟─d04690ca-390b-462b-8257-e9ebe01b3fd0
# ╟─e8201e38-73dc-45d4-b93d-20f063d2ac3a
# ╠═aa45acea-e6f7-49fb-a5b6-3ea0f2c1530c
# ╠═1e2848d6-252b-4db5-9864-1f1886da9998
# ╠═e2c3c7a8-0349-4f4c-87a7-3fbe8fc2fe46
# ╠═0e11d98d-79c6-444d-8df1-357e8218f233
# ╟─800262e6-5c78-4229-82d5-40de619d3013
# ╠═a6db5c0e-ea51-4071-b3fc-146ef90270aa
# ╠═ecf006d7-335c-48de-b629-b561dc08b334
# ╠═22fba19f-d337-48cc-a311-e44d03a0c050
# ╠═c41dc6d8-07ff-42fb-a54a-d1a1cf1cc223
# ╠═70e85333-26b7-40b8-8fe1-7fb470b5f7b3
# ╠═28f657cd-8ae7-4572-94fa-813797b59b25
