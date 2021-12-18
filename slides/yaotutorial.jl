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
using Pkg; ; Pkg.activate()

# ╔═╡ 2a145cba-26b0-43bd-9ab2-13818d246eae
using Revise, PlutoUI, Viznet, Compose

# ╔═╡ a8b06352-5265-4b89-bd82-b31f3cdac391
using Latexify

# ╔═╡ 57a3616e-49af-40b7-a000-a4ecc81af84e
using BitBasis

# ╔═╡ 342685aa-5159-11ec-13fd-fb8954106bca
using Yao, YaoPlots, Plots

# ╔═╡ 675b3398-01d7-4949-bb0e-7cdf9b805c69
using SymEngine: Basic  # import the symbolic data type

# ╔═╡ 0f8d63ac-f677-4889-b033-2a93f62be700
using YaoExtensions

# ╔═╡ 6f131f36-0b4e-4570-8527-620297fae48e
using YaoToEinsum, OMEinsumContractionOrders, OMEinsum

# ╔═╡ e11c26c0-e534-45fc-bb1c-c0f2ce4016db
SPACE = html"&nbsp; &nbsp; &nbsp; &nbsp;"

# ╔═╡ 0f099c85-f039-477e-a70d-a3801cbb2656
md"""
## Goals
* Introduce some open source packages for quantum simulation in Julia
* Full amplitude and tensor network based quantum simulation

## Contents
* Introduce Yao.jl (Roger Luo and Jinguo),
* From classical adder to reversible adder,
* What is new in Quantum, the Deutsch-Jozsa algorithm,
* Introduce OMEinsum (Andreas Peter and Jinguo)
* Towards faster simulation, a tensor network based quantum circuit simulation (basic).
"""

# ╔═╡ 48854a73-4896-4542-9ad4-15ae87418f1d
md"## Background knowledge"

# ╔═╡ c1d40103-1710-4221-b414-0958c13fb95f
md"""
I wish you are familiar with the following notations (~ 3rd year undergraduate)
* ``|\psi\rangle`` is a quantum state, or "ket",
* ``H`` is a Hamiltonian, it determines the dynamics of the quantum state as ``|\psi(t)\rangle = e^{-iHt}|\psi(0)\rangle``,
* ``\mathcal{O}`` is an observable, known as a Hermitian operator. One can measure this observable. The expectation value is the observed result is $\langle\mathcal{O}\rangle = \langle\psi|\mathcal{O}|\psi\rangle$
* ``X, Y`` and ``Z`` are Pauli operators, they anti-commute with each other: ``\{\sigma_i, \sigma_j\} = 2\delta_{ij}I`` for $\sigma_i\in \{X,Y,Z\}$.
"""

# ╔═╡ 4c8a4cea-63f4-49f0-82b9-dfa608be46bf
md"## If not..."

# ╔═╡ eb75810a-a746-4e4a-889f-c87c8d1e153f
html"""
<p>
<img src="https://static.docsity.com/documents_first_pages/notas/2012/04/17/dfe728fbf6999820eecb2c8fc0e773b6.png" width=300 style="vertical-align:sub;"><img src="https://images-na.ssl-images-amazon.com/images/I/51X+dIBIeZL._SY344_BO1,204,203,200_.jpg" width=290 style="vertical-align:super;">
</p>
"""

# ╔═╡ b83f7675-68fa-44b5-8681-c85984eeb877
md"# Quantum simulation with Yao"

# ╔═╡ 8017edf3-05b2-4cee-8a4e-90331a681037
let
	img = html"""<img src=https://docs.julialang.org/en/v1/assets/logo.svg width=100>"""
	md"""
## Why Julia?
$(img) is a **high-level**, **high-performance**, **dynamic programming language**. While it is a general-purpose language and can be used to write any application, many of its features are well suited for **numerical analysis** and **computational science**.

![](https://6lli539m39y3hpkelqsm3c2fg-wpengine.netdna-ssl.com/wp-content/uploads/2021/01/Julia-Adoption.png)
"""
end

# ╔═╡ 52d14933-99d4-4ee6-9eff-7be4a5722334
md"""## Why Yao?
Yao is an open source framework that aims to empower quantum information research with software tools. It is designed with following in mind:

* quantum algorithm design;
* quantum software 2.0;
* $(highlight{quantum computation education.})

by Roger Luo and Jinguo Liu, funded by Lei Wang and Pan Zhang.

$(html"<img src=https://camo.githubusercontent.com/477280de44a6d4408d3a3255d3d82a615a27ac1c5120063ef7d6b2f6640befb8/68747470733a2f2f79616f7175616e74756d2e6f72672f6173736574732f6c6f676f2e706e67 width=150>")
[https://yaoquantum.org/](https://yaoquantum.org/) (arXiv:1912.10877)


The most popular quantum simulator in Julia,
one of **Top 50** most popular Julia packages. It provides

* State of the art performance, with GPU support,
![](https://github.com/Roger-luo/quantum-benchmarks/raw/master/images/pcircuit.png)
* Matrix representation for operators,
* Built-in automatic differentiation engine,
* Generic data type (symbolic computation, tropical algera),
"""

# ╔═╡ b7d90c25-2f66-4ace-8d52-841ea376b3f9
md"## How classical gates work"

# ╔═╡ 2c6d024a-3187-4976-af27-393af8826a2d
md"This is an adder"

# ╔═╡ 165b369a-3d0a-4ae7-99f1-de0297f93707
md"![](http://www.worldofindie.co.uk/wp-content/uploads/2018/02/4-bit-adder-subtractor-circuit.png)"

# ╔═╡ 800bc5ce-1521-4c24-aaef-4911910d34cf
md"Adder in the configuration space"

# ╔═╡ c841f2e2-b907-4f74-be57-968ca339bec4
md"
FULL_ADDER = 
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
md"Univeral gate NAND and NOR"

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
md"NOT(1) = NAND(1, 1)"

# ╔═╡ c5694b66-d023-42c1-ae62-b7218ba8ebe7
md"AND(1, 2) = NOT(NAND(1, 2))"

# ╔═╡ c5812ca6-c4ca-4211-9d2d-df498fd7a2da
md"OR(1, 2) = NAND(NOT(1), NOT(2))"

# ╔═╡ 6b4e2611-5a5a-4646-aa8b-fec0d84d240c
md"## Reversible Gates"

# ╔═╡ 29d7eb5f-2d26-41af-80bc-1a14f0acce67
md"XOR(1, 2) is reversible"

# ╔═╡ 4591b3a0-ff3f-46c8-8c8d-ef7c6ba76cc5
md"NAND gate -> Toffoli Gate"

# ╔═╡ 657e1cfb-83ae-4f37-a237-15ff571e64b5
md"
```julia
001 → 001
011 → 011
101 → 101
111 → 110
000 → 000
010 → 010
100 → 100
110 → 111
```
"

# ╔═╡ 07f36b6e-a98d-4bee-9091-b1a717dab0e8
not_gate(n, i) = put(n, i=>X)

# ╔═╡ 7440555d-685f-4447-a44a-0463b37bd43c
xor_gate(n, i, j) = control(n, i, j=>X)

# ╔═╡ b2bf15bb-1350-4304-a932-c87b09558115
or_gate(n, i, j, k) = chain(kron(n, i=>X, j=>X, k=>X), control(n, (i, j), k=>X), kron(n, i=>X, j=>X))

# ╔═╡ 0e63c1eb-fea6-411a-bae6-0bc90dee6bc7
vizcircuit(or_gate(3, 1, 2, 3); show_ending_bar=false)

# ╔═╡ f311befa-ba09-42e6-ac9b-b59450162ebd
md"AND Gate"

# ╔═╡ 816ebc7b-43ec-47e2-a24a-86ac70dd6afe
and_gate(n, i, j, k) = control(n, (i, j), k=>X)

# ╔═╡ 893baeda-d5b0-43c1-9018-8d18c06486ca
vizcircuit(and_gate(3, 1, 2, 3); show_ending_bar=false)

# ╔═╡ 0c267c28-967b-4268-8372-00d7b48b7b8e
md"XOR gate"

# ╔═╡ f7b8dc5a-c589-41db-b550-d7b4c928dfef
full_adder() = xor_gate()

# ╔═╡ fed37413-c341-4dcf-a637-404d2f186b2b
vizcircuit(xor_gate(2, 1, 2); show_ending_bar=false)

# ╔═╡ 622fdaac-87dd-4c22-85b9-470510566480
# `s`, `cout` and ancillas (`x` and `y`) are initialized to 0
# flush data in `cin`
full_adder(n, a, b, cin, cout, x, y) = chain(and_gate(n, cin, a, x), xor_gate(n, a, cin), and_gate(n, b, cin, y), or_gate(n, x, y, cout), xor_gate(n, b, cin))

# ╔═╡ 26e5d619-7119-49bc-8907-17ae0db424f5
vizcircuit(full_adder(6, 1:6...); scale=0.7)

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
end

# ╔═╡ 617ba29e-a9f9-4ac4-83b6-232b2aff68d7
vizcircuit(add_circuit(2)[1]; scale=0.5)

# ╔═╡ e2e4ec60-2737-4560-89b1-1e14a35044e8
function calculate_binaryadd(n::Int, x, y)
	c, outputs = add_circuit(n)
	res = product_state(bit_literal([takebit(x, i) for i=1:n]..., [takebit(y, i) for i=1:n]..., zeros(Int, 3*n+1)...)) |> c |> measure!
	Int(readbit(res, outputs...))
end

# ╔═╡ aea490af-0bdd-4930-9ad2-7d9a13e08c46
calculate_binaryadd(2, 2, 2)

# ╔═╡ c0078012-3d81-4584-a050-9a58802d08a9
gatecount(add_circuit(4)[1])

# ╔═╡ ea458fa2-1f9f-46e1-88da-942034d0fa73
let
	n = 4
	circuit, out = add_circuit(n)
	vizcircuit(circuit; scale=0.3, starting_texts=[["a$i" for i=1:n]..., ["b$i" for i=1:n]..., ["0" for i=1:3n+1]...], ending_texts=[["a$i" for i=1:n]..., ["b$i" for i=1:n]..., [i ∈ out ? "c$(i-8)" : "" for i=1:3n+1]...])
end

# ╔═╡ 0840e7ea-aa63-44df-a788-ad18ac842006
md"## Balanced or constant?"

# ╔═╡ e5004508-2b9a-4acd-91ba-ad21fc7d1b16
md"""
``f: \{0, 1\}^n \rightarrow \{0, 1\}``

* *Balanced*: $f$ has an equal probability to output vaue 0 and 1.
* *Constant*: $f$ is a constant function.
"""

# ╔═╡ 2da83f06-e639-4f75-a162-cb561a8207ca
constantf(n) = chain(control(n+1, 1, n+1=>X), put(n+1, n+1=>X), control(n+1, 1, n+1=>X))

# ╔═╡ 408f534d-5553-4a95-8d14-cabf36f8aff5
vizcircuit(constantf(4), starting_texts=[["a$i" for i=1:4]..., "0"], ending_texts=[["a$i" for i=1:4]..., "f"])

# ╔═╡ 774def48-8f9f-425a-aa47-ffe396692247
md"## Quantum bits and quantum gates"

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
deutsch_jozsa(circuit::AbstractBlock{N}, inputs, output) where N = chain(put(N+1, N+1=>X), repeat(N+1, H, [inputs..., N+1]), subroutine(N+1, circuit, 1:N), control(N+1, output, N+1=>X),  repeat(N+1, H, inputs), Yao.Measure(N+1, locs=inputs))

# ╔═╡ d8e9ae38-1569-497a-8506-a3e059ffc7ab
let
	c = deutsch_jozsa(constantf(4), 1:4, 5)
	vizcircuit(c)
end

# ╔═╡ 195ccd7e-f1ea-4704-9392-ffe596a756f9
let
	addc, outputs = add_circuit(4)
	c = deutsch_jozsa(addc, 1:4, outputs[1])
	vizcircuit(c; scale=0.3)
end

# ╔═╡ f1814dfa-588d-49b8-a68e-bda185f3712f
let
	c = deutsch_jozsa(constantf(4), 1:4, 5)
	zero_state(nqubits(c)) |> c
	c[end].results
end

# ╔═╡ 023a5424-ab33-433a-bcdf-2315735aa00d
let
	addc, outputs = add_circuit(4)
	c = deutsch_jozsa(addc, 1:4, outputs[1])
	zero_state(nqubits(c)) |> c
	c[end].results
end

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

# ╔═╡ 1e2848d6-252b-4db5-9864-1f1886da9998
vizcircuit(decompose_toffoli(control(3, (1,2), 3=>X)); w_depth=0.7, scale=0.6)

# ╔═╡ e2c3c7a8-0349-4f4c-87a7-3fbe8fc2fe46
vizcircuit(decompose_toffoli(add_circuit(4)[1]); scale=0.3)

# ╔═╡ 0e11d98d-79c6-444d-8df1-357e8218f233
gatecount(decompose_toffoli(add_circuit(4)[1]))

# ╔═╡ ee4186d7-1f93-4b72-a86c-076b136e2eda
md"## Not gate"

# ╔═╡ 0f655562-0c18-43a8-83f0-4d11b7eeb068
vizcircuit_with_inputs(g, inputs) = vizcircuit(g; show_ending_bar=false, starting_texts=inputs, starting_offset=-0.3, ending_texts=measure!(product_state(bit_literal(inputs...)) |> g), ending_offset=0.3);

# ╔═╡ 9f402b0e-f6ae-4a2c-960e-5b756a1f6739
vizcircuit_with_inputs(X, [1])

# ╔═╡ 8dac807a-445e-477a-8dca-b717277039e0
vizcircuit_with_inputs(X, [0])

# ╔═╡ e188093d-1ef1-4e70-966b-77cc0761a801
Latexify.LaTeXString("X = " * latexify(mat(Basic, X); env=:raw).s)

# ╔═╡ 02e018b3-2f2e-42b7-b5ab-e06bc72daad4
cnot_gate = control(2, 1, 2=>X);

# ╔═╡ 13fdcd91-4360-4245-8cb6-dbc0b2238976
vizcircuit(cnot_gate; show_ending_bar=false, starting_texts=[1, 1], starting_offset=-0.3, ending_texts=[0, 1], ending_offset=0.3)

# ╔═╡ 67ad78c7-930a-48b1-b8d3-d9e549ee7379
latexify(mat(Basic, cnot_gate))

# ╔═╡ 38194214-3bf1-4229-9fe8-37282b30a5ad
latexify(mat(Basic, control(3, (3,2), 1=>X)))

# ╔═╡ cee8db91-7e84-4728-aef6-ea861c62ff96
vizcircuit(control(3, (1,2), 3=>X); show_ending_bar=false, starting_texts=[1,1,0], starting_offset=-0.3, ending_texts=[1,1,1], ending_offset=0.3)

# ╔═╡ b6a356a8-de90-405f-ad5e-6ae9131f871e
md"""## The linear algebra representation of Quantum machanics
Linear algebra to quantum dynamics is what regular number to abstract algebra.
"""

# ╔═╡ 400f697d-fbdb-4067-881f-e004f50bbc0f
md"A quantum state is represented as a **normalized** vector. e.g. the following vector represents a quantum system with 5 states ($s=1,2,\ldots,5$)"

# ╔═╡ b3daf0a3-4d1e-4474-bb99-948c718059ce
#   normalized  random-gaussian   complex-valued  vector
ψ0 = normalize!(   randn(          ComplexF64,      5))

# ╔═╡ 6a719cf5-d8a4-4369-95e7-c1dfbc5cebf8
md"""
Its elements ``\langle s|\psi\rangle =\sqrt{p_s}e^{i \phi_s}``. $p$ is the probability given by Born's rule.
"""

# ╔═╡ 681cdc3b-348f-45c4-9031-6298e81d995d
md"""
Hamiltonian is a **Hermitian** matrix
"""

# ╔═╡ b3813f2d-4e9d-4789-a653-d28fd1b98eda
rand_hami = rand_hermitian(5)

# ╔═╡ 6294a05c-46b1-47e4-9f3c-e741ab191a09
obs = rand_hermitian(5)

# ╔═╡ bf983162-66e3-40e0-afce-7f5b7728c490
@bind evolve_time Slider(0:0.03:2; show_value=true)

# ╔═╡ d85a2028-9b6f-4d6c-ba61-fecef3fd4876
ψt = exp(-im*rand_hami*evolve_time) * ψ0

# ╔═╡ 0139dbd2-486e-4b58-b656-6b5e06864cd1
mean_obs = ψt' * obs * ψt

# ╔═╡ fe47198a-5f63-4171-90d2-ce1e9e5ec0a2
md"## Qubits"

# ╔═╡ 30ef3d84-f85c-45f3-9256-f413637ccb2f
md"A random state"

# ╔═╡ def1dae7-cd10-4d75-9a7b-cbe02b430609
md"""
⚫ $SPACE  ⚪ $SPACE ⚫ $SPACE ⚫ $SPACE ⚪ $SPACE ⚫

``~s_6``  ``~~~~~~~s_5`` ``~~~~~~~s_4`` ``~~~~~~~s_3`` ``~~~~~~~s_2`` ``~~~~~~~s_1`` 
"""

# ╔═╡ edd1809e-9087-4301-bbea-31420699d777
reg_rand = rand_state(6)

# ╔═╡ 21c0c2d7-40e6-4ec5-88b7-3af0b8774c41
statevec(reg_rand)

# ╔═╡ 5f28cbbe-c8b9-424b-9b40-7f5babc4fc8f
md"⚪ $SPACE  ⚪ $SPACE ⚪ $SPACE ⚪ $SPACE ⚪ $SPACE ⚪"

# ╔═╡ 8d1c023b-4ee6-4201-99f6-0245cee70cda
reg_zero = zero_state(6)

# ╔═╡ bab73814-016d-4040-a994-553070d03ae7
md"⚪ $SPACE  ⚪ $SPACE ⚫ $SPACE ⚪ $SPACE ⚪ $SPACE ⚪"

# ╔═╡ 6f64057b-130d-44c1-8957-e4c7c56f7cb4
reg_product = product_state(bit"001000")

# ╔═╡ 30100576-9fa2-4825-bfa6-35ae6cb920bb
md"## Operators"

# ╔═╡ bfb9fcb5-2555-490a-9dcc-48cc91f8ab3c
md"""
#### Primitive Operators
"""

# ╔═╡ 82b6f24c-4eab-4420-848e-66026d5aa8ee
typeof(X)

# ╔═╡ 2e17662d-24ed-40c9-93b0-4cef526c3a75
@bind selected_gate Select([X, Y, Z, I2, ConstGate.P0, ConstGate.P1, ConstGate.Pu, ConstGate.Pd, ConstGate.T, ConstGate.S, SWAP, Yao.Measure(1)])

# ╔═╡ 21b11fac-5efb-4fd5-a1b0-e684d215a46c
vizcircuit(selected_gate)

# ╔═╡ a2129e65-7e73-4b42-9924-e88d60893ed2
!(selected_gate isa Yao.Measure) && mat(Basic, selected_gate)

# ╔═╡ 8b1606b5-78d8-4020-8597-1cf9f2364a2b
mat(rot(X, Basic(:θ)))

# ╔═╡ 8be8d92a-d882-4e77-a3c3-9846babeda2b
mat(rot(SWAP, Basic(:θ)))

# ╔═╡ 597bff8b-f7fa-4fc8-bc59-f06c934f211e
md"""
#### Composite Operators
"""

# ╔═╡ 421abcef-a1be-4297-a159-adcd92332d6b
md"
`put(n, (i,j...,)=>G)`
```math
I_1\otimes I_2 \otimes \ldots G_{i,j,\ldots} \ldots \otimes I_n
```
"

# ╔═╡ 3388350c-fd8b-4cae-bf34-71446bd747ca
g_put1 = put(2, 2=>selected_gate)

# ╔═╡ 99847765-fd49-4c9d-914b-5d47eac89490
vizcircuit(g_put1)

# ╔═╡ fbd1ea6e-16da-4289-b0e3-d8e8a5a55541
!(selected_gate isa Yao.Measure) && mat(Basic, g_put1)

# ╔═╡ 5dafbbbd-d6f7-423c-8a46-41b1596f83a3
g_put2 = put(3, (3, 1)=>SWAP)

# ╔═╡ 984db7f2-30e0-421d-9f27-874553b6fae2
mat(Basic, g_put2)

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
g_ctrl1 = control(2, 2, 1=>selected_gate)

# ╔═╡ 3aaa7e49-9507-4384-ac7e-f4d486127811
vizcircuit(g_ctrl1)

# ╔═╡ 65a01ff1-8458-4198-990c-46814a16b714
!(selected_gate isa Yao.Measure) && mat(Basic, g_ctrl1)

# ╔═╡ 47da74c7-0f24-4dd5-b6f7-88c22a81de5a
g_ctrl2 = control(2, -2, 1=>selected_gate)

# ╔═╡ e06b23aa-7fae-4495-9b69-43cea3682a28
vizcircuit(g_ctrl2)

# ╔═╡ 3e154273-d97f-4c71-a500-48913766c836
!(selected_gate isa Yao.Measure) && mat(Basic, g_ctrl2)

# ╔═╡ 915185d6-3b36-437e-bac0-22cd5503bba2
mat(Basic, control(3, (3,2), 1=>X))

# ╔═╡ 0a6d5a41-f42b-42e7-adf7-93a40f67c9f3
md"`chain(G1, G2, ..., Gn)` or `Gn * ... * G_2 * G_1`
```math
G_n G_{n-1} \ldots G_1
```
"

# ╔═╡ dc46cbbe-0eb5-41f9-ad5f-f8a7fba5ff9c
mat(Basic, chain(X, Y))

# ╔═╡ 0a7cab72-a50b-4d58-877a-4a9d939b198f
g_chain2 = chain(put(2, 1=>selected_gate), control(2, 1, 2=>X))

# ╔═╡ 6b2d8bb9-5a4f-4699-ac92-bc8a4003918f
vizcircuit(g_chain2)

# ╔═╡ 26435992-fd7c-4797-8fc9-9c5f52fc6388
!(selected_gate isa Yao.Measure) && mat(Basic, g_chain2)

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
heisenberg(n) = sum([sum([kron(n, i=>G, i+1=>G) for G in [X, Y, Z]]) for i=1:n-1])

# ╔═╡ 3fa73e8a-fe46-475f-8346-e54f52c144f4
md"Hermitian: ``\mathcal{O} = \mathcal{O}^\dagger``"

# ╔═╡ a4b4761b-0686-42bd-9f04-d6083b6064e4
heisenberg(3)

# ╔═╡ 31c1c4fc-7476-4066-8462-bf18f8d69966
ishermitian(heisenberg(3))

# ╔═╡ 81a55867-4158-4371-a12c-b08e3d64c0cb
md"Unitary: ``\mathcal{O}^\dagger\mathcal{O} = I``"

# ╔═╡ 7f15cbb0-09be-46ee-852c-e43dde9bc4f5
isunitary(heisenberg(3))

# ╔═╡ cc4e7d7b-0596-4bc0-b23e-f6861fcd5260
md"Reflexive: ``\mathcal{O}^2 = I``"

# ╔═╡ 8b03f2f1-f38c-4d5e-b778-0c3b4aaf910d
isreflexive(heisenberg(3))

# ╔═╡ f639e1b4-e257-4218-8727-292eadedd7bf
mat(heisenberg(20))

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

# ╔═╡ be280ec7-c42c-413b-b4de-274600190c17
md"## Time evolution"

# ╔═╡ c5e36e92-9717-45bc-91e0-568ae206d57c
te = time_evolve(heisenberg(5), 0.3)

# ╔═╡ dce2db9d-c9ff-451d-ba91-0143e7a116cc
apply!(rand_state(5), te)

# ╔═╡ 17960ae3-8ba7-4d70-90a1-a161a988a4de
md"## Measure"

# ╔═╡ 64e71a26-991a-43d8-9d4b-73920ac1b532
reg_measure = zero_state(5)

# ╔═╡ fe69d633-1da8-4db5-9e30-eb19c4ec392f
reg_measure |> dispatch!(variational_circuit(5, 2), :random)

# ╔═╡ f02a5014-45ed-4bfd-817f-4544ca15f0a0
measure_results = measure(reg_measure; nshots=1000);

# ╔═╡ 23a5b68f-bc93-43a6-9072-7d929f5f62dc
bar(0:31, [count(b->Int(b)==i, measure_results) for i=0:31], label="counting")

# ╔═╡ 840e56f4-2fcf-4136-bdb7-7feba88f25da
md"## Tensor network based quantum simulation"

# ╔═╡ 09137c47-1bd0-481b-8b60-ca2abff64afc
md"""## Vector is not the only choice
```math
v := \left(\begin{matrix}a_1\\ b_1\end{matrix}\right) \otimes \left(\begin{matrix}a_2\\ b_2\end{matrix}\right) \otimes \ldots \otimes \left(\begin{matrix}a_n\\ b_n\end{matrix}\right)
```
"""

# ╔═╡ 62393bfb-5597-40cd-a5e4-c4b10a61d1cf
md"""
#### References:
* Simulating quantum computation by contracting tensor networks, Igor L. Markov, Yaoyun Shi (arXiv:quant-ph/0511069)
* Solving the sampling problem of the Sycamore quantum supremacy circuits, Feng Pan, Keyang Chen, Pan Zhang (arXiv:2111.03011)
"""

# ╔═╡ fdf13c95-6d95-4220-9759-5fea44eea274
optcode = let
	c = decompose_toffoli(add_circuit(4)[1])
	code, xs = yao2einsum(c; initial_state=Dict([i=>rand(0:1) for i=1:nqubits(c)]))
	optcode = optimize_code(code, uniformsize(code, 2), TreeSA())
end

# ╔═╡ 669bed15-1c57-401b-bcfb-48b5a1520a58
timespace_complexity(optcode, uniformsize(optcode, 2))

# ╔═╡ d04690ca-390b-462b-8257-e9ebe01b3fd0
md"""
## Discussion
* Github issue (Yao.jl)
* Julia slack (channel: yao-dev)
"""

# ╔═╡ Cell order:
# ╠═4a96f5c9-37b4-4a8a-a6bd-8a4b4440eb49
# ╠═2a145cba-26b0-43bd-9ab2-13818d246eae
# ╠═e11c26c0-e534-45fc-bb1c-c0f2ce4016db
# ╠═a8b06352-5265-4b89-bd82-b31f3cdac391
# ╟─0f099c85-f039-477e-a70d-a3801cbb2656
# ╟─48854a73-4896-4542-9ad4-15ae87418f1d
# ╟─c1d40103-1710-4221-b414-0958c13fb95f
# ╟─4c8a4cea-63f4-49f0-82b9-dfa608be46bf
# ╟─eb75810a-a746-4e4a-889f-c87c8d1e153f
# ╟─b83f7675-68fa-44b5-8681-c85984eeb877
# ╟─8017edf3-05b2-4cee-8a4e-90331a681037
# ╠═52d14933-99d4-4ee6-9eff-7be4a5722334
# ╟─b7d90c25-2f66-4ace-8d52-841ea376b3f9
# ╟─2c6d024a-3187-4976-af27-393af8826a2d
# ╟─165b369a-3d0a-4ae7-99f1-de0297f93707
# ╟─800bc5ce-1521-4c24-aaef-4911910d34cf
# ╟─c841f2e2-b907-4f74-be57-968ca339bec4
# ╟─e1b6b7a9-2d26-4f43-a0bd-54f7ad22a5b3
# ╟─42f3e654-1836-491a-8119-b03b93822f45
# ╟─cc0d4b1c-07b6-42ed-b1e2-80581a02ee6f
# ╟─c5694b66-d023-42c1-ae62-b7218ba8ebe7
# ╟─c5812ca6-c4ca-4211-9d2d-df498fd7a2da
# ╟─6b4e2611-5a5a-4646-aa8b-fec0d84d240c
# ╟─29d7eb5f-2d26-41af-80bc-1a14f0acce67
# ╟─4591b3a0-ff3f-46c8-8c8d-ef7c6ba76cc5
# ╟─657e1cfb-83ae-4f37-a237-15ff571e64b5
# ╠═07f36b6e-a98d-4bee-9091-b1a717dab0e8
# ╠═7440555d-685f-4447-a44a-0463b37bd43c
# ╠═b2bf15bb-1350-4304-a932-c87b09558115
# ╠═0e63c1eb-fea6-411a-bae6-0bc90dee6bc7
# ╟─f311befa-ba09-42e6-ac9b-b59450162ebd
# ╠═816ebc7b-43ec-47e2-a24a-86ac70dd6afe
# ╠═893baeda-d5b0-43c1-9018-8d18c06486ca
# ╟─0c267c28-967b-4268-8372-00d7b48b7b8e
# ╠═f7b8dc5a-c589-41db-b550-d7b4c928dfef
# ╠═fed37413-c341-4dcf-a637-404d2f186b2b
# ╠═622fdaac-87dd-4c22-85b9-470510566480
# ╠═26e5d619-7119-49bc-8907-17ae0db424f5
# ╠═a7d85b82-a705-4f3b-a371-06a87071335d
# ╠═617ba29e-a9f9-4ac4-83b6-232b2aff68d7
# ╠═e2e4ec60-2737-4560-89b1-1e14a35044e8
# ╠═aea490af-0bdd-4930-9ad2-7d9a13e08c46
# ╠═c0078012-3d81-4584-a050-9a58802d08a9
# ╠═ea458fa2-1f9f-46e1-88da-942034d0fa73
# ╟─0840e7ea-aa63-44df-a788-ad18ac842006
# ╟─e5004508-2b9a-4acd-91ba-ad21fc7d1b16
# ╟─2da83f06-e639-4f75-a162-cb561a8207ca
# ╟─408f534d-5553-4a95-8d14-cabf36f8aff5
# ╟─774def48-8f9f-425a-aa47-ffe396692247
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
# ╠═023a5424-ab33-433a-bcdf-2315735aa00d
# ╟─e8201e38-73dc-45d4-b93d-20f063d2ac3a
# ╠═aa45acea-e6f7-49fb-a5b6-3ea0f2c1530c
# ╠═1e2848d6-252b-4db5-9864-1f1886da9998
# ╠═e2c3c7a8-0349-4f4c-87a7-3fbe8fc2fe46
# ╠═0e11d98d-79c6-444d-8df1-357e8218f233
# ╟─ee4186d7-1f93-4b72-a86c-076b136e2eda
# ╠═57a3616e-49af-40b7-a000-a4ecc81af84e
# ╠═0f655562-0c18-43a8-83f0-4d11b7eeb068
# ╟─9f402b0e-f6ae-4a2c-960e-5b756a1f6739
# ╟─8dac807a-445e-477a-8dca-b717277039e0
# ╠═e188093d-1ef1-4e70-966b-77cc0761a801
# ╟─02e018b3-2f2e-42b7-b5ab-e06bc72daad4
# ╟─13fdcd91-4360-4245-8cb6-dbc0b2238976
# ╠═67ad78c7-930a-48b1-b8d3-d9e549ee7379
# ╟─38194214-3bf1-4229-9fe8-37282b30a5ad
# ╠═cee8db91-7e84-4728-aef6-ea861c62ff96
# ╟─b6a356a8-de90-405f-ad5e-6ae9131f871e
# ╟─400f697d-fbdb-4067-881f-e004f50bbc0f
# ╠═b3daf0a3-4d1e-4474-bb99-948c718059ce
# ╟─6a719cf5-d8a4-4369-95e7-c1dfbc5cebf8
# ╟─681cdc3b-348f-45c4-9031-6298e81d995d
# ╠═b3813f2d-4e9d-4789-a653-d28fd1b98eda
# ╠═6294a05c-46b1-47e4-9f3c-e741ab191a09
# ╠═bf983162-66e3-40e0-afce-7f5b7728c490
# ╠═d85a2028-9b6f-4d6c-ba61-fecef3fd4876
# ╠═0139dbd2-486e-4b58-b656-6b5e06864cd1
# ╟─fe47198a-5f63-4171-90d2-ce1e9e5ec0a2
# ╠═342685aa-5159-11ec-13fd-fb8954106bca
# ╟─30ef3d84-f85c-45f3-9256-f413637ccb2f
# ╟─def1dae7-cd10-4d75-9a7b-cbe02b430609
# ╠═edd1809e-9087-4301-bbea-31420699d777
# ╠═21c0c2d7-40e6-4ec5-88b7-3af0b8774c41
# ╟─5f28cbbe-c8b9-424b-9b40-7f5babc4fc8f
# ╠═8d1c023b-4ee6-4201-99f6-0245cee70cda
# ╟─bab73814-016d-4040-a994-553070d03ae7
# ╠═6f64057b-130d-44c1-8957-e4c7c56f7cb4
# ╟─30100576-9fa2-4825-bfa6-35ae6cb920bb
# ╟─bfb9fcb5-2555-490a-9dcc-48cc91f8ab3c
# ╠═82b6f24c-4eab-4420-848e-66026d5aa8ee
# ╠═2e17662d-24ed-40c9-93b0-4cef526c3a75
# ╠═21b11fac-5efb-4fd5-a1b0-e684d215a46c
# ╠═a2129e65-7e73-4b42-9924-e88d60893ed2
# ╠═8b1606b5-78d8-4020-8597-1cf9f2364a2b
# ╠═8be8d92a-d882-4e77-a3c3-9846babeda2b
# ╟─597bff8b-f7fa-4fc8-bc59-f06c934f211e
# ╟─421abcef-a1be-4297-a159-adcd92332d6b
# ╠═3388350c-fd8b-4cae-bf34-71446bd747ca
# ╠═675b3398-01d7-4949-bb0e-7cdf9b805c69
# ╠═99847765-fd49-4c9d-914b-5d47eac89490
# ╠═fbd1ea6e-16da-4289-b0e3-d8e8a5a55541
# ╠═5dafbbbd-d6f7-423c-8a46-41b1596f83a3
# ╠═984db7f2-30e0-421d-9f27-874553b6fae2
# ╟─f35b9789-788b-450d-af8e-472412efce77
# ╠═34a0b380-ac20-4d88-b498-969745a0df7c
# ╠═3aaa7e49-9507-4384-ac7e-f4d486127811
# ╠═65a01ff1-8458-4198-990c-46814a16b714
# ╠═47da74c7-0f24-4dd5-b6f7-88c22a81de5a
# ╠═e06b23aa-7fae-4495-9b69-43cea3682a28
# ╠═3e154273-d97f-4c71-a500-48913766c836
# ╠═915185d6-3b36-437e-bac0-22cd5503bba2
# ╟─0a6d5a41-f42b-42e7-adf7-93a40f67c9f3
# ╠═dc46cbbe-0eb5-41f9-ad5f-f8a7fba5ff9c
# ╠═0a7cab72-a50b-4d58-877a-4a9d939b198f
# ╠═6b2d8bb9-5a4f-4699-ac92-bc8a4003918f
# ╠═26435992-fd7c-4797-8fc9-9c5f52fc6388
# ╟─3a89d593-ae7e-454c-bb18-0baa4cc06a25
# ╠═4559f70b-5a9b-478d-ba9a-3a35f6fc9716
# ╟─6cac6c4b-a0d7-4cae-bc2d-1cf9e0d775a2
# ╠═1b35ea7b-c831-4fce-8109-3cd35a80214a
# ╟─3fa73e8a-fe46-475f-8346-e54f52c144f4
# ╠═a4b4761b-0686-42bd-9f04-d6083b6064e4
# ╠═31c1c4fc-7476-4066-8462-bf18f8d69966
# ╟─81a55867-4158-4371-a12c-b08e3d64c0cb
# ╠═7f15cbb0-09be-46ee-852c-e43dde9bc4f5
# ╟─cc4e7d7b-0596-4bc0-b23e-f6861fcd5260
# ╠═8b03f2f1-f38c-4d5e-b778-0c3b4aaf910d
# ╠═f639e1b4-e257-4218-8727-292eadedd7bf
# ╟─800262e6-5c78-4229-82d5-40de619d3013
# ╠═a6db5c0e-ea51-4071-b3fc-146ef90270aa
# ╠═ecf006d7-335c-48de-b629-b561dc08b334
# ╠═22fba19f-d337-48cc-a311-e44d03a0c050
# ╠═c41dc6d8-07ff-42fb-a54a-d1a1cf1cc223
# ╠═70e85333-26b7-40b8-8fe1-7fb470b5f7b3
# ╠═28f657cd-8ae7-4572-94fa-813797b59b25
# ╟─be280ec7-c42c-413b-b4de-274600190c17
# ╠═c5e36e92-9717-45bc-91e0-568ae206d57c
# ╠═dce2db9d-c9ff-451d-ba91-0143e7a116cc
# ╟─17960ae3-8ba7-4d70-90a1-a161a988a4de
# ╠═0f8d63ac-f677-4889-b033-2a93f62be700
# ╠═64e71a26-991a-43d8-9d4b-73920ac1b532
# ╠═fe69d633-1da8-4db5-9e30-eb19c4ec392f
# ╠═f02a5014-45ed-4bfd-817f-4544ca15f0a0
# ╠═23a5b68f-bc93-43a6-9072-7d929f5f62dc
# ╟─840e56f4-2fcf-4136-bdb7-7feba88f25da
# ╟─09137c47-1bd0-481b-8b60-ca2abff64afc
# ╟─62393bfb-5597-40cd-a5e4-c4b10a61d1cf
# ╠═6f131f36-0b4e-4570-8527-620297fae48e
# ╠═fdf13c95-6d95-4220-9759-5fea44eea274
# ╠═669bed15-1c57-401b-bcfb-48b5a1520a58
# ╟─d04690ca-390b-462b-8257-e9ebe01b3fd0
