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

# ╔═╡ 342685aa-5159-11ec-13fd-fb8954106bca
using Yao, YaoPlots, Plots

# ╔═╡ 675b3398-01d7-4949-bb0e-7cdf9b805c69
using SymEngine: Basic  # import the symbolic data type

# ╔═╡ 0f8d63ac-f677-4889-b033-2a93f62be700
using YaoExtensions

# ╔═╡ 6f131f36-0b4e-4570-8527-620297fae48e
using YaoToEinsum

# ╔═╡ e11c26c0-e534-45fc-bb1c-c0f2ce4016db
SPACE = html"&nbsp; &nbsp; &nbsp; &nbsp;"

# ╔═╡ 0f099c85-f039-477e-a70d-a3801cbb2656
md"""
## Goals
* quantum simulation basics,
* high level understanding of tensor network based quantum circuit simulation,
* Julia (a programming language) quantum toolbox
    * Yao
    * OMEinsum
"""

# ╔═╡ 48854a73-4896-4542-9ad4-15ae87418f1d
md"## Background knowledge"

# ╔═╡ c1d40103-1710-4221-b414-0958c13fb95f
md"""
I wish you are familiar with the following notations
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

# ╔═╡ b6a356a8-de90-405f-ad5e-6ae9131f871e
md"""## The linear algebra representation of Quantum machanics
Linear algebra to quantum dynamics is what regular number to abstract algebra.
"""

# ╔═╡ 400f697d-fbdb-4067-881f-e004f50bbc0f
md"A quantum state is represented as a **normalized** vector. e.g. the following vector represents a quantum system with 5 states ($s=1,2,\ldots,5$)"

# ╔═╡ b3daf0a3-4d1e-4474-bb99-948c718059ce
ψ0 = normalize!(randn(ComplexF64, 5))

# ╔═╡ 6a719cf5-d8a4-4369-95e7-c1dfbc5cebf8
md"""
Its elements ``\langle s|\psi\rangle =\sqrt{p_s}e^{i \phi_s}``.
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
Yao is a ... by Roger Luo and Jinguo Liu, funded by Lei Wang and Pan Zhang.

$(html"<img src=https://camo.githubusercontent.com/477280de44a6d4408d3a3255d3d82a615a27ac1c5120063ef7d6b2f6640befb8/68747470733a2f2f79616f7175616e74756d2e6f72672f6173736574732f6c6f676f2e706e67 width=200>")
[https://yaoquantum.org/](https://yaoquantum.org/)

The most popular quantum simulator in Julia,
one of **Top 50** most popular Julia packages. It provides

* State of the art performance, with GPU support,
![](https://github.com/Roger-luo/quantum-benchmarks/raw/master/images/pcircuit.png)
* Matrix representation for operators,
* Built-in automatic differentiation engine,
* Generic data type (symbolic computation, tropical algera),
"""

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
vizcircuit(selected_gate; show_ending=false)

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
vizcircuit(g_put1; show_ending=false)

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
vizcircuit(g_ctrl1; show_ending=false)

# ╔═╡ 65a01ff1-8458-4198-990c-46814a16b714
!(selected_gate isa Yao.Measure) && mat(Basic, g_ctrl1)

# ╔═╡ 47da74c7-0f24-4dd5-b6f7-88c22a81de5a
g_ctrl2 = control(2, -2, 1=>selected_gate)

# ╔═╡ e06b23aa-7fae-4495-9b69-43cea3682a28
vizcircuit(g_ctrl2; show_ending=false)

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
vizcircuit(g_chain2; show_ending=false)

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
	YaoPlots.CircuitStyles.linecolor[] = "#ffffff"
	YaoPlots.CircuitStyles.textcolor[] = "#ffffff"
	YaoPlots.CircuitStyles.gate_bgcolor[] = "white"
	vizcircuit(qft_circuit(5))
end |> PNG("qft_white.png")

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
# ╟─0f099c85-f039-477e-a70d-a3801cbb2656
# ╟─48854a73-4896-4542-9ad4-15ae87418f1d
# ╟─c1d40103-1710-4221-b414-0958c13fb95f
# ╟─4c8a4cea-63f4-49f0-82b9-dfa608be46bf
# ╟─eb75810a-a746-4e4a-889f-c87c8d1e153f
# ╟─b6a356a8-de90-405f-ad5e-6ae9131f871e
# ╠═342685aa-5159-11ec-13fd-fb8954106bca
# ╟─400f697d-fbdb-4067-881f-e004f50bbc0f
# ╠═b3daf0a3-4d1e-4474-bb99-948c718059ce
# ╟─6a719cf5-d8a4-4369-95e7-c1dfbc5cebf8
# ╟─681cdc3b-348f-45c4-9031-6298e81d995d
# ╠═b3813f2d-4e9d-4789-a653-d28fd1b98eda
# ╠═6294a05c-46b1-47e4-9f3c-e741ab191a09
# ╠═bf983162-66e3-40e0-afce-7f5b7728c490
# ╠═d85a2028-9b6f-4d6c-ba61-fecef3fd4876
# ╠═0139dbd2-486e-4b58-b656-6b5e06864cd1
# ╟─b83f7675-68fa-44b5-8681-c85984eeb877
# ╟─8017edf3-05b2-4cee-8a4e-90331a681037
# ╟─52d14933-99d4-4ee6-9eff-7be4a5722334
# ╟─fe47198a-5f63-4171-90d2-ce1e9e5ec0a2
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
# ╟─d04690ca-390b-462b-8257-e9ebe01b3fd0
