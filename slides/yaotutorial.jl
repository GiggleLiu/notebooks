### A Pluto.jl notebook ###
# v0.17.3

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

# ╔═╡ 2a145cba-26b0-43bd-9ab2-13818d246eae
using Revise, PlutoUI, Viznet, Compose, Latexify

# ╔═╡ 57a3616e-49af-40b7-a000-a4ecc81af84e
using BitBasis

# ╔═╡ 40563837-b1e9-4df9-9cc7-8841d0068973
using Yao, YaoPlots

# ╔═╡ 8faad9e1-8b29-4721-be29-fab2ce3c0e4c
using SymEngine: Basic  # the symbolics type

# ╔═╡ 5fbddde1-3594-43c4-9f56-c0ae258d926f
using Yao.ConstGate: P0

# ╔═╡ 6f131f36-0b4e-4570-8527-620297fae48e
using YaoToEinsum, OMEinsumContractionOrders, OMEinsum

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

# ╔═╡ 620f31d6-45ec-470f-bf21-0eee32214666
md"A gate set can form any logical expression."

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

# ╔═╡ 3381001d-1120-4b88-ac01-5ca861f0a9be
md"Note: ancillas must be initialized to a known value, e.g. 0"

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

# ╔═╡ ea458fa2-1f9f-46e1-88da-942034d0fa73
let
	n = 4
	circuit, out = add_circuit(n)
	vizcircuit(circuit; scale=0.3, starting_texts=[["a$i" for i=1:n]..., ["b$i" for i=1:n]..., ["0" for i=1:3n+1]...], ending_texts=[["a$i" for i=1:n]..., ["b$i" for i=1:n]..., [i+2n ∈ out ? "c$(i)" : "?" for i=1:3n+1]...])
end

# ╔═╡ 22cd8a36-c76f-4a8c-a8a2-e1924136012a
md"## Adder with uncomputing"

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
	vizcircuit(circuit; scale=0.15, starting_texts=[["a$i" for i=1:n]..., ["b$i" for i=1:n]..., ["0" for i=1:4n+1]...], ending_texts=[["a$i" for i=1:n]..., ["b$i" for i=1:n]..., [i+2n ∈ out ? "c$(i)" : "0" for i=1:4n+1]...])
end

# ╔═╡ b69f7d3f-cc92-4e9d-813d-55643a2a3b30
md"""This is $(highlight("compute-copy-uncompute")), it can bring polynomial overhead in time/space. """

# ╔═╡ cb2c496c-e9dc-4666-b41b-d97cca377047
md"## Universal reversible gate

* {Toffoli}
* {Fredkin}
"

# ╔═╡ 176f3757-9f03-4990-840e-2cebccd6abaa
md"The combination of these gate can form any permutaion matrix."

# ╔═╡ 78cb38db-b0b6-49f7-bcaf-1a2df58b8150
md"CNOT = "

# ╔═╡ 334db854-7805-4548-a414-ad7d215fe387
vizcircuit(and_gate(3, 1, 2, 3); starting_texts=["1", "a", "b"])

# ╔═╡ 40956c0e-3833-4397-a4a8-c46e1890ef39
md"NOT = "

# ╔═╡ 8bebf527-d197-48b0-9778-f71b0a5dad77
vizcircuit(and_gate(3, 1, 2, 3); starting_texts=["1", "1", "a"])

# ╔═╡ ca094f81-bef3-40eb-9279-086fe7eb506a
md"""## Irreversible v.s. Reversible
* Reversible computing may have polynomial time/space overhead.
* Reversible computing can be much more $(highlight("energy efficient")) (arxiv: 1803.02789), because erasing 1 bit information costs a least ``k_bT`` energy (known as the Landauer's principle).
"""

# ╔═╡ 0840e7ea-aa63-44df-a788-ad18ac842006
md"## Balanced or constant?"

# ╔═╡ e5004508-2b9a-4acd-91ba-ad21fc7d1b16
md"""
``f: \{0, 1\}^n \rightarrow \{0, 1\}``

<<<<<<< HEAD
# ╔═╡ 861d41c6-c5d0-45cc-9557-5dfce41e68a3
YaoToEinsum.yao2einsum(decompose_toffoli(add_circuit(4)[1]))

# ╔═╡ ee4186d7-1f93-4b72-a86c-076b136e2eda
md"## Not gate"
=======
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
    * Toffoli gate
    * CNOT gate
* One hot representation of a classical state, which has a probabilistic interpretation. A reversible gate is a permutation matrix.
* Univeral gate sets
* The problem of distinguishing balanced or constant functions.

### TODO
* A tutorial of Yao,
* Build a quantum adder with Yao,
* Use a quantum algorithm to solve the "balanced or constant" problem,
* Tensor network based simulation of quantum circuits,
"""
>>>>>>> 1636c0ca50e7f532eec60625f2ff877ec7aae9c8

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
Hamiltonian is a **Hermitian** matrix, its time evolve operator is unitary
"""

# ╔═╡ 45d19896-ae93-4f9a-b1fe-bb143c7ff580
md"
```math
|\psi(t)\rangle = e^{-iHt} |\psi(0)\rangle
```
"

# ╔═╡ 4fd68654-bb38-499f-8e93-1872a3ded4db
md"Permutation matrix is a special unitary matrix."

# ╔═╡ f4f18b91-86f6-460c-b5de-4c52931a4098
md"## Quantum universal gate set
* \{Toffoli, H\}

* \{CNOT, H, S, T\}
"

# ╔═╡ 2b231cca-a450-4773-867a-65c2d367ed1d
md"The combination of these gate can form any unitary matrix."

# ╔═╡ bfb9fcb5-2555-490a-9dcc-48cc91f8ab3c
md"""
## Primitive Operators
"""

# ╔═╡ 2e17662d-24ed-40c9-93b0-4cef526c3a75
@bind selected_gate Select([X, Y, Z, I2, ConstGate.P0, ConstGate.P1, ConstGate.Pu, ConstGate.Pd, ConstGate.T, ConstGate.S, SWAP, ConstGate.CNOT, ConstGate.Toffoli, Yao.Measure(1)])

# ╔═╡ db743480-4c67-465d-89a8-f0d2e6d0d152
md"Note: for PDF version, this selection box might not work. For the original notebook, check the link on the top of this file."

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
vizcircuit(g_ctrl1)   # Pluto notebook/VSCode

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

# ╔═╡ c2e846cb-9f80-4ea1-8295-eb9dadd7b878
md"Measuring on computational basis ``\{P_1^{(n)}, P_1^{(n-1)}, \ldots P_1^{(1)}\}``"

# ╔═╡ 8a94551d-c82a-4db0-aaf0-3c5285ab0695
let
	reg = rand_state(2)
	measure!(reg), statevec(reg)
end

# ╔═╡ 9eb46218-30c9-4109-9b32-9382fdee1081
let
	reg = rand_state(2)
	measure!(obs, reg), statevec(reg)
end

# ╔═╡ 957b8318-a4a6-45aa-839d-4f301ccf3c80
mat(obs)

# ╔═╡ 847051a1-48e9-4852-8f6d-9a5f604e9142
md"## Deutsch–Jozsa algorithm"

# ╔═╡ 85228d6e-731d-465d-a000-672e7fff8aff
md"![](https://upload.wikimedia.org/wikipedia/commons/thumb/b/b5/Deutsch-Jozsa-algorithm-quantum-circuit.png/600px-Deutsch-Jozsa-algorithm-quantum-circuit.png)"

# ╔═╡ f66bd812-cfd0-4ff2-be21-71c38cb8ab79
md"Case 1, ``f(x)`` is a constant.
"

# ╔═╡ 6d6b55d4-4ddb-4ef9-b2b1-154021c01273
md"Case 2, ``f(x)`` is balanced."

# ╔═╡ 5467dafc-265a-4848-a898-364124a6e9b4
mat(chain(H, H))    # hadamard gate is reflexive

# ╔═╡ 1e01ef7f-63d0-4b0d-8455-7b524206f15f
statevec(zero_state(1) |> X |> H)   # `|>` is the pipline operator

# ╔═╡ 1ca545b9-6f88-4251-aa8a-d5ad0236c77e
statevec(zero_state(1) |> X |> H |> X)

# ╔═╡ 528bfa77-d6da-4700-a72c-5afffdbe6139
zero_state(4) |> repeat(4, H, 1:4) |> statevec

# ╔═╡ 1571a8ae-aae5-4af8-9b9c-6fe9e832a39e
md"""
```math
\begin{align}
|\psi\rangle &= \sum_x|x\rangle \otimes \left(X^{f(x)}\frac{|0\rangle - |1\rangle}{\sqrt{2}}\right)\\
 &= \sum_{x\in\{x|f(x)=0\}}|x\rangle \otimes \frac{|0\rangle - |1\rangle}{\sqrt{2}} - \sum_{x\in\{x|f(x)=1\}}|x\rangle \otimes \frac{|0\rangle - |1\rangle}{\sqrt{2}}\\
&= \underbrace{\sum_{x\in\{x|f(x)=0\}}|x\rangle - \sum_{x\in\{x|f(x)=1\}}|x\rangle}_{\text{has zero overlap with $(H |0\rangle)^{\otimes n}$}} \otimes \frac{|0\rangle - |1\rangle}{\sqrt{2}}
\end{align}
```
"""

# ╔═╡ c60ea349-6b96-4f70-8369-a92379338bbd
deutsch_jozsa(circuit::AbstractBlock{N}, inputs, output) where N = chain(
	put(N+1, N+1=>X),                   # initialize the last qubit to state `1`
	repeat(N+1, H, [inputs..., N+1]),   # hadamard gates
	subroutine(N+1, circuit, 1:N),      # compute
	control(N+1, output, N+1=>X),       # CNOT used for copy
	subroutine(N+1, circuit', 1:N),     # uncompute
	repeat(N+1, H, inputs),             # hadamard gates
	Yao.Measure(N+1, locs=inputs)       # measure
)

# ╔═╡ 6c370435-4bca-42f6-975e-b633b5611444
md"# Quantum adder?
The same as reversible circuits."

# ╔═╡ 79d21733-6d95-4d95-a337-261fb483f4f0
# `s`, `cout` and ancillas (`x` and `y`) are initialized to 0
# flush data in `cin`
q_full_adder(n, a, b, cin, cout, x, y) = chain(and_gate(n, cin, a, x), xor_gate(n, a, cin), and_gate(n, b, cin, y), or_gate(n, x, y, cout), xor_gate(n, b, cin));

# ╔═╡ ba3ac294-5f74-4169-966f-f8d93b48253b
q_constantf(n) = chain(control(n+1, 1, n+1=>X), put(n+1, n+1=>X), control(n+1, 1, n+1=>X));

# ╔═╡ 82119d4b-a08a-4c0e-b4ee-6a018e88731e
vizcircuit(q_constantf(4), starting_texts=[["a$i" for i=1:4]..., "0"], ending_texts=[["a$i" for i=1:4]..., "f"])

# ╔═╡ 129135d6-9623-4e4f-812a-fcfe15d1e5f5
vizcircuit(q_full_adder(6, 1:6...); scale=0.7, starting_texts=["a", "b", "cᵢ", "0", "0", "0"], ending_texts=["a", "b", "cᵢ'", "cₒ", "?", "?"])

# ╔═╡ dc345279-84ca-4a4e-b37d-74a644e0a83a
function q_add_circuit(n::Int)
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
		blk = q_full_adder(nbit, as[i], bs[i], cin, cout, xs[i], ys[i])
		push!(c, blk)
		cs[i] = cin
		cin = cout
	end
	cs[end] = cin
	return c, cs
end;

# ╔═╡ 74dc8276-35c3-43a1-a92b-b7116fef6bb1
function q_adder_with_uncomputing(n)
	add, outputs = q_add_circuit(n)
	n1 = nqubits(add)
	chain(subroutine(n1+n, add, 1:n1), [control(n1+n, outputs[i], n1+i=>X) for i=1:n]..., subroutine(n1+n, add', 1:n1)), collect(n1+1:n1+n)
end;

# ╔═╡ 63917a68-597e-40e5-a479-624fc80d7cc6
vizcircuit(q_adder_with_uncomputing(4)[1]; scale=0.3)

# ╔═╡ 8d3c9ef3-1dec-4369-aa8c-ecff882ace6b
md"## Use Deutsch-Jozsa algorithm to solve the balance-constant problem."

# ╔═╡ d8e9ae38-1569-497a-8506-a3e059ffc7ab
let
	c = deutsch_jozsa(q_constantf(4), 1:4, 5)
	vizcircuit(c; scale=0.5)
end

# ╔═╡ 11f74d18-1fb4-405f-8923-195a861029f4
md"The balanced function"

# ╔═╡ 195ccd7e-f1ea-4704-9392-ffe596a756f9
let
	addc, outputs = q_add_circuit(2)
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
statevec(zero_state(4) |> repeat(4, H, 1:4))

# ╔═╡ 84d4d1fd-0921-48ea-b292-fa2ccce58e7d
statevec(zero_state(4) |> repeat(4, H, 1:4) |> put(4, 1=>Z))

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

<<<<<<< HEAD
# ╔═╡ e5acf299-875b-44cd-80dc-7dbb20c1c4f3

=======
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

# ╔═╡ ac82ea9c-be86-4501-95db-04f7eef26bea
cluster_code = let
	c = mps_like_circuit(10)
	code, xs = yao2einsum(c; initial_state=Dict([i=>zero_state(1) for i=1:nqubits(c)]), final_state=Dict([i=>zero_state(1) for i=1:nqubits(c)]))
	optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
end

# ╔═╡ b931b25c-508e-4e64-bc3a-ff993643e165
timespace_complexity(cluster_code, uniformsize(cluster_code, 2))

# ╔═╡ fdf13c95-6d95-4220-9759-5fea44eea274
optcode = let
	c = add_circuit(4)[1]
	code, xs = yao2einsum(c; initial_state=Dict([i=>rand(0:1) for i=1:nqubits(c)]), final_state=Dict([i=>rand(0:1) for i=1:nqubits(c)]))
	optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
end

# ╔═╡ 669bed15-1c57-401b-bcfb-48b5a1520a58
timespace_complexity(optcode, uniformsize(optcode, 2))
>>>>>>> 1636c0ca50e7f532eec60625f2ff877ec7aae9c8

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
* jinguoliu@g.harvard.edu
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BitBasis = "50ba71b6-fa0f-514d-ae9a-0916efc90dcf"
Compose = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
OMEinsum = "ebe7aa44-baf0-506c-a96f-8464559b3922"
OMEinsumContractionOrders = "6f22d1fd-8eed-4bb7-9776-e7d684900715"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Revise = "295af30f-e4ad-537b-8983-00126c2a3abe"
SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
Viznet = "52a3aca4-6234-47fd-b74a-806bdf78ede9"
Yao = "5872b779-8223-5990-8dd0-5abbb0748c8c"
YaoPlots = "32cfe2d9-419e-45f2-8191-2267705d8dbc"
YaoToEinsum = "9b173c7b-dc24-4dc5-a0e1-adab2f7b6ba9"

[compat]
BitBasis = "~0.7.4"
Compose = "~0.9.2"
Latexify = "~0.15.9"
OMEinsum = "~0.6.5"
OMEinsumContractionOrders = "~0.6.2"
PlutoUI = "~0.7.23"
Revise = "~3.2.1"
SymEngine = "~0.8.7"
Viznet = "~0.3.3"
Yao = "~0.6.4"
YaoPlots = "~0.7.0"
YaoToEinsum = "~0.1.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0"
manifest_format = "2.0"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "485ee0867925449198280d4af84bdb46a2a404d0"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.0.1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "abb72771fd8895a7ebd83d5632dc4b989b022b5b"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "f87e559f87a45bece9c9ed97458d3afe98b1ebb9"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.1.0"

[[deps.ArrayInterface]]
deps = ["Compat", "IfElse", "LinearAlgebra", "Requires", "SparseArrays", "Static"]
git-tree-sha1 = "265b06e2b1f6a216e0e8f183d28e4d354eab3220"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "3.2.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.BFloat16s]]
deps = ["LinearAlgebra", "Printf", "Random", "Test"]
git-tree-sha1 = "a598ecb0d717092b5539dbbe890c98bac842b072"
uuid = "ab4f0b2a-ad5b-11e8-123f-65d77653426b"
version = "0.2.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BatchedRoutines]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "8ee75390ba4bbfaf9aa48c121857b0da9a914265"
uuid = "a9ab73d0-e05c-5df1-8fde-d6a4645b8d8e"
version = "0.2.1"

[[deps.BetterExp]]
git-tree-sha1 = "dd3448f3d5b2664db7eceeec5f744535ce6e759b"
uuid = "7cffe744-45fd-4178-b173-cf893948b8b7"
version = "0.1.0"

[[deps.BitBasis]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "68ce92be119ad7ff44ebbb9ffc0f7a70b1e34c45"
uuid = "50ba71b6-fa0f-514d-ae9a-0916efc90dcf"
version = "0.7.4"

[[deps.CEnum]]
git-tree-sha1 = "215a9aa4a1f23fbd05b92769fdd62559488d70e9"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.1"

[[deps.CUDA]]
deps = ["AbstractFFTs", "Adapt", "BFloat16s", "CEnum", "CompilerSupportLibraries_jll", "ExprTools", "GPUArrays", "GPUCompiler", "LLVM", "LazyArtifacts", "Libdl", "LinearAlgebra", "Logging", "Printf", "Random", "Random123", "RandomNumbers", "Reexport", "Requires", "SparseArrays", "SpecialFunctions", "TimerOutputs"]
git-tree-sha1 = "2c8329f16addffd09e6ca84c556e2185a4933c64"
uuid = "052768ef-5323-5732-b1bb-66c8b64840ba"
version = "3.5.0"

[[deps.CacheServers]]
deps = ["Distributed", "Test"]
git-tree-sha1 = "b584b04f236d3677b4334fab095796a128445bf8"
uuid = "a921213e-d44a-5460-ac04-5d720a99ba71"
version = "0.2.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "4c26b4e9e91ca528ea212927326ece5918a04b47"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.2"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "9aa8a5ebb6b5bf469a7e0e2b5202cf6f8c291104"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.0.6"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "32a2b8af383f11cbb65803883837a149d10dfe8a"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.10.12"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Compose]]
deps = ["Base64", "Colors", "DataStructures", "Dates", "IterTools", "JSON", "LinearAlgebra", "Measures", "Printf", "Random", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "c6461fc7c35a4bb8d00905df7adafcff1fe3a6bc"
uuid = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
version = "0.9.2"

[[deps.DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Dierckx]]
deps = ["Dierckx_jll"]
git-tree-sha1 = "5fefbe52e9a6e55b8f87cb89352d469bd3a3a090"
uuid = "39dd38d3-220a-591b-8e3c-4c3a8c710a94"
version = "0.5.1"

[[deps.Dierckx_jll]]
deps = ["CompilerSupportLibraries_jll", "Libdl", "Pkg"]
git-tree-sha1 = "a580560f526f6fc6973e8bad2b036514a4e3b013"
uuid = "cd4c43a9-7502-52ba-aa6d-59fb2a88580b"
version = "0.0.1+0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.ExponentialUtilities]]
deps = ["ArrayInterface", "LinearAlgebra", "Printf", "Requires", "SparseArrays"]
git-tree-sha1 = "1b873816d2cfc8c0fcb1edcb08e67fdf630a70b7"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.10.2"

[[deps.ExprTools]]
git-tree-sha1 = "b7e3d17636b348f005f11040025ae8c6f645fe92"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.6"

[[deps.Expronicon]]
deps = ["MLStyle", "Pkg", "TOML"]
git-tree-sha1 = "eb43b420c63be5df51549ae86c43ed72f20ebcb9"
uuid = "6b7a57c9-7cc1-4fdf-b7f5-e857abae3636"
version = "0.6.13"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.GMP_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "781609d7-10c4-51f6-84f2-b8444358ff6d"

[[deps.GPUArrays]]
deps = ["Adapt", "LinearAlgebra", "Printf", "Random", "Serialization", "Statistics"]
git-tree-sha1 = "7772508f17f1d482fe0df72cabc5b55bec06bbe0"
uuid = "0c68f7d7-f131-5f86-a1c3-88cf8149b2d7"
version = "8.1.2"

[[deps.GPUCompiler]]
deps = ["ExprTools", "InteractiveUtils", "LLVM", "Libdl", "Logging", "TimerOutputs", "UUIDs"]
git-tree-sha1 = "2cac236070c2c4b36de54ae9146b55ee2c34ac7a"
uuid = "61eb1bfa-7361-4325-ad38-22787b887f55"
version = "0.13.10"

[[deps.GraphPlot]]
deps = ["ArnoldiMethod", "ColorTypes", "Colors", "Compose", "DelimitedFiles", "Graphs", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "5e51d9d9134ebcfc556b82428521fe92f709e512"
uuid = "a2cc645c-3eea-5389-862e-a155d0052231"
version = "0.5.0"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "92243c07e786ea3458532e199eb3feee0e7e08eb"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.4.1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "e273807f38074f033d94207a201e6e827d8417db"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.8.21"

[[deps.LLVM]]
deps = ["CEnum", "LLVMExtra_jll", "Libdl", "Printf", "Unicode"]
git-tree-sha1 = "7cc22e69995e2329cc047a879395b2b74647ab5f"
uuid = "929cbde3-209d-540e-8aea-75f648917ca0"
version = "4.7.0"

[[deps.LLVMExtra_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c5fc4bef251ecd37685bea1c4068a9cfa41e8b9a"
uuid = "dad2f222-ce93-54a1-a47d-0025e8a3acab"
version = "0.0.13+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a8f4f279b6fa3c3c4f1adadd78a621b13a506bce"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.9"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LegibleLambdas]]
deps = ["MacroTools"]
git-tree-sha1 = "7946db4829eb8de47c399f92c19790f9cc0bbd07"
uuid = "f1f30506-32fe-5131-bd72-7c197988f9e5"
version = "0.3.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "e5718a00af0ab9756305a0392832c8952c7426c1"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.6"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "491a883c4fef1103077a7f648961adbf9c8dd933"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.1.2"

[[deps.LuxurySparse]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "f6eb35c3a10571c1d62748250197c66ed4c42607"
uuid = "d05aeea4-b7d4-55ac-b691-9e7fabb07ba2"
version = "0.6.8"

[[deps.MLStyle]]
git-tree-sha1 = "594e189325f66e23a8818e5beb11c43bb0141bcd"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.10"

[[deps.MPC_jll]]
deps = ["Artifacts", "GMP_jll", "JLLWrappers", "Libdl", "MPFR_jll", "Pkg"]
git-tree-sha1 = "9618bed470dcb869f944f4fe4a9e76c4c8bf9a11"
uuid = "2ce0c516-f11f-5db3-98ad-e0e1048fbd70"
version = "1.2.1+0"

[[deps.MPFR_jll]]
deps = ["Artifacts", "GMP_jll", "Libdl"]
uuid = "3a97d323-0669-5f0c-9066-3539efd106a3"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.Multigraphs]]
deps = ["Graphs", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "055a7c49a626e17a8c99bcaaf472d0de60848929"
uuid = "7ebac608-6c66-46e6-9856-b5f43e107bac"
version = "0.3.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OMEinsum]]
deps = ["AbstractTrees", "BatchedRoutines", "CUDA", "ChainRulesCore", "Combinatorics", "LinearAlgebra", "MacroTools", "Requires", "Test", "TupleTools"]
git-tree-sha1 = "c172922074434ef8dda952da7178208ad832637a"
uuid = "ebe7aa44-baf0-506c-a96f-8464559b3922"
version = "0.6.5"

[[deps.OMEinsumContractionOrders]]
deps = ["BetterExp", "OMEinsum", "Requires", "SparseArrays", "Suppressor"]
git-tree-sha1 = "c3f853756b1c1b52d2959ebe3921366b724c367d"
uuid = "6f22d1fd-8eed-4bb7-9776-e7d684900715"
version = "0.6.2"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "d7fa6237da8004be601e19bd6666083056649918"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.3"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "5152abbdab6488d5eec6a01029ca6697dff4ec8f"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.23"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Random123]]
deps = ["Libdl", "Random", "RandomNumbers"]
git-tree-sha1 = "0e8b146557ad1c6deb1367655e052276690e71a3"
uuid = "74087812-796a-5b5d-8853-05524746bad3"
version = "1.4.2"

[[deps.RandomNumbers]]
deps = ["Random", "Requires"]
git-tree-sha1 = "043da614cc7e95c703498a491e2c21f58a2b8111"
uuid = "e6cf234a-135c-5ec9-84dd-332b85af5143"
version = "1.5.3"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "8f82019e525f4d5c669692772a6f4b0a58b06a6a"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.2.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "e55f4c73ec827f96cd52db0bc6916a3891c726b5"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "f0bccf98e16759818ffc5d97ac3ebf87eb950150"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.8.1"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "7f5a513baec6f122401abfc8e9c074fdac54f6c1"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.4.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
git-tree-sha1 = "0f2aa8e32d511f758a2ce49208181f7733a0936a"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.1.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "2bb0cb32026a66037360606510fca5984ccc6b75"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.13"

[[deps.Suppressor]]
git-tree-sha1 = "a819d77f31f83e5792a76081eee1ea6342ab8787"
uuid = "fd094767-a336-5f1f-9728-57cf17d0bbfb"
version = "0.2.0"

[[deps.SymEngine]]
deps = ["Compat", "Libdl", "LinearAlgebra", "RecipesBase", "SpecialFunctions", "SymEngine_jll"]
git-tree-sha1 = "6cf88a0b98c758a36e6e978a41e8a12f6f5cdacc"
uuid = "123dc426-2d89-5057-bbad-38513e3affd8"
version = "0.8.7"

[[deps.SymEngine_jll]]
deps = ["Artifacts", "GMP_jll", "JLLWrappers", "Libdl", "MPC_jll", "MPFR_jll", "Pkg"]
git-tree-sha1 = "3cd0f249ae20a0093f839738a2f2c1476d5581fe"
uuid = "3428059b-622b-5399-b16f-d347a77089a4"
version = "0.8.1+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "7cb456f358e8f9d102a8b25e8dfedf58fa5689bc"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.13"

[[deps.TupleTools]]
git-tree-sha1 = "3c712976c47707ff893cf6ba4354aa14db1d8938"
uuid = "9d95972d-f1c8-5527-a6e0-b4b365fa01f6"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Viznet]]
deps = ["Compose", "Dierckx"]
git-tree-sha1 = "7a022ae6ac8b153d47617ed8c196ce60645689f1"
uuid = "52a3aca4-6234-47fd-b74a-806bdf78ede9"
version = "0.3.3"

[[deps.Yao]]
deps = ["BitBasis", "Reexport", "YaoArrayRegister", "YaoBase", "YaoBlocks", "YaoSym"]
git-tree-sha1 = "21e49c3b1f3ec891fd728664feae600c8601b013"
uuid = "5872b779-8223-5990-8dd0-5abbb0748c8c"
version = "0.6.4"

[[deps.YaoAPI]]
git-tree-sha1 = "dc4edfcda2e59fd2624f84941da040a4e30220e3"
uuid = "0843a435-28de-4971-9e8b-a9641b2983a8"
version = "0.1.0"

[[deps.YaoArrayRegister]]
deps = ["Adapt", "BitBasis", "LinearAlgebra", "LuxurySparse", "Random", "StaticArrays", "StatsBase", "TupleTools", "YaoBase"]
git-tree-sha1 = "0a29506643daf3d39f38b337d42ac287ffb9b535"
uuid = "e600142f-9330-5003-8abb-0ebd767abc51"
version = "0.7.10"

[[deps.YaoBase]]
deps = ["BitBasis", "LegibleLambdas", "LinearAlgebra", "LuxurySparse", "MLStyle", "Random", "Reexport", "SparseArrays", "TupleTools", "YaoAPI"]
git-tree-sha1 = "bd91eaf91b5a5c5a8091b5ade97810dc63695851"
uuid = "a8f54c17-34bc-5a9d-b050-f522fe3f755f"
version = "0.14.5"

[[deps.YaoBlocks]]
deps = ["BitBasis", "CacheServers", "ChainRulesCore", "ExponentialUtilities", "InteractiveUtils", "LegibleLambdas", "LinearAlgebra", "LuxurySparse", "MLStyle", "Random", "SimpleTraits", "SparseArrays", "StaticArrays", "StatsBase", "TupleTools", "YaoArrayRegister", "YaoBase"]
git-tree-sha1 = "e3f0a62b1934dd75ae751b0d73ef99f9ffa0ceed"
uuid = "418bc28f-b43b-5e0b-a6e7-61bbc1a2c1df"
version = "0.11.10"

[[deps.YaoHIR]]
deps = ["Expronicon", "MLStyle", "YaoLocations"]
git-tree-sha1 = "938fb1436c702dddbbfb2bf173f28e8ebbc48ee0"
uuid = "6769671a-fce8-4286-b3f7-6099e1b1298a"
version = "0.2.0"

[[deps.YaoLocations]]
git-tree-sha1 = "c90c42c8668c9096deb0c861822f0f8f80cbdc68"
uuid = "66df03fb-d475-48f7-b449-3d9064bf085b"
version = "0.1.6"

[[deps.YaoPlots]]
deps = ["BitBasis", "Colors", "Compose", "GraphPlot", "Graphs", "Multigraphs", "Viznet", "YaoBlocks", "ZXCalculus"]
git-tree-sha1 = "bae26fcc84ae8082755dcdd087c943fc06b5aebd"
uuid = "32cfe2d9-419e-45f2-8191-2267705d8dbc"
version = "0.7.0"

[[deps.YaoSym]]
deps = ["BitBasis", "LinearAlgebra", "LuxurySparse", "Requires", "SparseArrays", "YaoArrayRegister", "YaoBase", "YaoBlocks"]
git-tree-sha1 = "e34838fa98d02d4c969ba9f92783a12a336e2f88"
uuid = "3b27209a-d3d6-11e9-3c0f-41eb92b2cb9d"
version = "0.4.7"

[[deps.YaoToEinsum]]
deps = ["LinearAlgebra", "OMEinsum", "Yao"]
git-tree-sha1 = "2b22e59f0f73b0fa978085269e6d944ee6864808"
uuid = "9b173c7b-dc24-4dc5-a0e1-adab2f7b6ba9"
version = "0.1.2"

[[deps.ZXCalculus]]
deps = ["Graphs", "LinearAlgebra", "MLStyle", "Multigraphs", "SparseArrays", "YaoHIR", "YaoLocations"]
git-tree-sha1 = "58e4f9a72618f2daf483f328fd82f0d10df8dc37"
uuid = "3525faa3-032d-4235-a8d4-8c2939a218dd"
version = "0.5.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─2a145cba-26b0-43bd-9ab2-13818d246eae
# ╟─e6306a69-bd6a-4c01-9c6e-1cb523668019
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
# ╟─620f31d6-45ec-470f-bf21-0eee32214666
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
<<<<<<< HEAD
# ╠═f7b8dc5a-c589-41db-b550-d7b4c928dfef
# ╠═fed37413-c341-4dcf-a637-404d2f186b2b
# ╠═622fdaac-87dd-4c22-85b9-470510566480
# ╠═26e5d619-7119-49bc-8907-17ae0db424f5
# ╠═a7d85b82-a705-4f3b-a371-06a87071335d
# ╠═617ba29e-a9f9-4ac4-83b6-232b2aff68d7
# ╠═e2e4ec60-2737-4560-89b1-1e14a35044e8
# ╠═aea490af-0bdd-4930-9ad2-7d9a13e08c46
# ╠═c0078012-3d81-4584-a050-9a58802d08a9
# ╠═aa45acea-e6f7-49fb-a5b6-3ea0f2c1530c
# ╠═1e2848d6-252b-4db5-9864-1f1886da9998
# ╠═008d228b-5d88-44c2-a4bc-c877329349bd
# ╠═e2c3c7a8-0349-4f4c-87a7-3fbe8fc2fe46
# ╠═0e11d98d-79c6-444d-8df1-357e8218f233
# ╠═861d41c6-c5d0-45cc-9557-5dfce41e68a3
# ╟─ee4186d7-1f93-4b72-a86c-076b136e2eda
# ╠═57a3616e-49af-40b7-a000-a4ecc81af84e
# ╠═0f655562-0c18-43a8-83f0-4d11b7eeb068
# ╟─9f402b0e-f6ae-4a2c-960e-5b756a1f6739
# ╟─8dac807a-445e-477a-8dca-b717277039e0
# ╠═e188093d-1ef1-4e70-966b-77cc0761a801
# ╟─02e018b3-2f2e-42b7-b5ab-e06bc72daad4
# ╟─13fdcd91-4360-4245-8cb6-dbc0b2238976
# ╠═67ad78c7-930a-48b1-b8d3-d9e549ee7379
=======
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
>>>>>>> 1636c0ca50e7f532eec60625f2ff877ec7aae9c8
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
# ╟─3381001d-1120-4b88-ac01-5ca861f0a9be
# ╟─a7d85b82-a705-4f3b-a371-06a87071335d
# ╟─86b398d0-844d-4d06-9ede-79d618502cce
# ╟─e2e4ec60-2737-4560-89b1-1e14a35044e8
# ╟─aea490af-0bdd-4930-9ad2-7d9a13e08c46
# ╟─7f84eed8-edd8-4f8e-a2b1-3ad862285934
# ╟─ea458fa2-1f9f-46e1-88da-942034d0fa73
# ╟─22cd8a36-c76f-4a8c-a8a2-e1924136012a
# ╟─051548ff-dd39-4d55-ac53-e8e2bacec68e
# ╟─e7b2e327-952e-4062-9193-19653eaee19c
# ╟─b69f7d3f-cc92-4e9d-813d-55643a2a3b30
# ╟─cb2c496c-e9dc-4666-b41b-d97cca377047
# ╟─176f3757-9f03-4990-840e-2cebccd6abaa
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
# ╟─45d19896-ae93-4f9a-b1fe-bb143c7ff580
# ╟─4fd68654-bb38-499f-8e93-1872a3ded4db
# ╟─f4f18b91-86f6-460c-b5de-4c52931a4098
# ╟─2b231cca-a450-4773-867a-65c2d367ed1d
# ╟─bfb9fcb5-2555-490a-9dcc-48cc91f8ab3c
# ╟─2e17662d-24ed-40c9-93b0-4cef526c3a75
# ╟─db743480-4c67-465d-89a8-f0d2e6d0d152
# ╠═21b11fac-5efb-4fd5-a1b0-e684d215a46c
# ╟─3fa73e8a-fe46-475f-8346-e54f52c144f4
# ╠═31c1c4fc-7476-4066-8462-bf18f8d69966
# ╟─81a55867-4158-4371-a12c-b08e3d64c0cb
# ╠═7f15cbb0-09be-46ee-852c-e43dde9bc4f5
# ╟─cc4e7d7b-0596-4bc0-b23e-f6861fcd5260
# ╠═8b03f2f1-f38c-4d5e-b778-0c3b4aaf910d
# ╠═8faad9e1-8b29-4721-be29-fab2ce3c0e4c
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
# ╟─c2e846cb-9f80-4ea1-8295-eb9dadd7b878
# ╠═8a94551d-c82a-4db0-aaf0-3c5285ab0695
# ╠═9eb46218-30c9-4109-9b32-9382fdee1081
# ╠═957b8318-a4a6-45aa-839d-4f301ccf3c80
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
# ╟─6c370435-4bca-42f6-975e-b633b5611444
# ╠═79d21733-6d95-4d95-a337-261fb483f4f0
# ╠═ba3ac294-5f74-4169-966f-f8d93b48253b
# ╠═82119d4b-a08a-4c0e-b4ee-6a018e88731e
# ╠═129135d6-9623-4e4f-812a-fcfe15d1e5f5
# ╠═dc345279-84ca-4a4e-b37d-74a644e0a83a
# ╠═74dc8276-35c3-43a1-a92b-b7116fef6bb1
# ╠═63917a68-597e-40e5-a479-624fc80d7cc6
# ╟─8d3c9ef3-1dec-4369-aa8c-ecff882ace6b
# ╠═d8e9ae38-1569-497a-8506-a3e059ffc7ab
# ╟─11f74d18-1fb4-405f-8923-195a861029f4
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
<<<<<<< HEAD
# ╠═e5acf299-875b-44cd-80dc-7dbb20c1c4f3
=======
# ╠═7a39257d-b864-460d-bccd-eed9933c4b7c
# ╠═19c37ef5-fe1d-4061-a5ad-f4b6bd9ecba4
# ╠═ac82ea9c-be86-4501-95db-04f7eef26bea
# ╠═b931b25c-508e-4e64-bc3a-ff993643e165
# ╠═fdf13c95-6d95-4220-9759-5fea44eea274
# ╠═669bed15-1c57-401b-bcfb-48b5a1520a58
>>>>>>> 1636c0ca50e7f532eec60625f2ff877ec7aae9c8
# ╟─d04690ca-390b-462b-8257-e9ebe01b3fd0
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
