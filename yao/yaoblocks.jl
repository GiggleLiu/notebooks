### A Pluto.jl notebook ###
# v0.19.8

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

# ╔═╡ 173125da-e761-11ec-25b9-fd328db9847e
using Yao, YaoPlots, PlutoUI, BenchmarkTools, KrylovKit

# ╔═╡ 7eedb3d5-a3b6-40e5-8d2d-94623faf8d92
TableOfContents()

# ╔═╡ 84705dc3-dd2c-419b-983e-e6c879dd79da
md"# YaoBlocks - Hypercubic Linear algebra"

# ╔═╡ 1c835e36-299d-41c7-940f-6e8d9fea87f9
html"""
<div align="center">
<a class="Header-link " href="https://github.com/QuantumBFS/Yao.jl" data-hotkey="g d" aria-label="Homepage " data-ga-click="Header, go to dashboard, icon:logo">
  <svg class="octicon octicon-mark-github v-align-middle" height="32" viewBox="0 0 16 16" version="1.1" width="32" aria-hidden="true"><path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"></path></svg>
</a>
<br>
<a href="https://raw.githubusercontent.com/GiggleLiu/notebooks/master/notebooks/yao/yaoblocks.jl" target="_blank"> download this notebook </a></div>
"""

# ╔═╡ caec900e-13fd-4340-8c18-5c1d4ee91a29
md"## Definition and interfaces
Unlike most quantum simulators, an Yao's block system is not only for a quantum gate, it implements the **linear algebra on a hypercubic**.
All blocks subtypes `AbstractBlock`, and are classified into two categories
1. `PrimitiveBlock`: a leaf block that does not contain sub-blocks.
2. `CompositeBlock`: a block composed of other blocks.
"

# ╔═╡ 6b5100df-e813-4ed8-8f7c-c4fa05015ea6
@with_terminal print_blocktree()

# ╔═╡ 2d7486af-c5e7-465a-b8d9-d8f7b00d4db4
qft3 = EasyBuild.qft_circuit(3)

# ╔═╡ c63152ea-9a68-4cff-bacf-978ef95ea9a8
vizcircuit(qft3)

# ╔═╡ 31619f1a-784b-48c4-9755-ff41d273f46f
md"### The `mat` interface"

# ╔═╡ 23a0f167-481e-4de4-b175-b96ee1923d7d
md"""The only "must have" interface is `mat`."""

# ╔═╡ 63e8fea9-01c7-4f78-a493-acff4f412820
mat(H)

# ╔═╡ 0afe664a-2ce5-4867-b6c2-8169d492c4bc
md"""
`PutBlock` is a composite block, its matrix representation is:
```math
\texttt{put(n, i=>G)} = I^{\otimes n-i} \otimes G \otimes I^{\otimes i-1}
```
"""

# ╔═╡ f2391ac6-8bc3-40c3-bb08-4604e4fd8dda
vizcircuit(put(3, 1=>H))

# ╔═╡ 221f1d4d-814a-4d51-9aa5-1a963acebf23
mat(put(3, 1=>H))

# ╔═╡ f274ba5b-6a4a-40bf-af55-1a0bf10119a4
md"""
```math
{\rm shift}(\theta) = \begin{pmatrix}
1 & 0\\
0 & e^{i\theta}
\end{pmatrix}
```
"""

# ╔═╡ d53ef0be-74c5-4f0e-820f-62ef4455f340
mat(shift(π/2))

# ╔═╡ 3a739d22-f847-4459-bf1f-f7a58817a0ae
md"Control block is a composite block, its matrix representation is"

# ╔═╡ 729a8205-657d-4b1b-92a3-750ffe4ececb
md"""
```math
\begin{align}
\texttt{control(n, i, j=>G)} = &I^{\otimes n-i} P_0 \otimes I^{\otimes i-j-1} \otimes I\otimes I^{\otimes j-1}
+\\
& I^{\otimes n-i} P_1 \otimes I^{\otimes i-j-1} \otimes G\otimes I^{\otimes j-1}
\end{align}
```
"""

# ╔═╡ de0ea7f9-cd37-4e36-9bf3-dcc18658ad74
md"""``P_0=|0\rangle\langle 0|``"""

# ╔═╡ 0b4c480d-8c72-4437-b18b-7c87a3abdf32
mat(ConstGate.P0)

# ╔═╡ 05063b40-e87f-49d8-ac15-66371b8e8cce
md"""``P_1=|1\rangle\langle 1|``"""

# ╔═╡ ef0f8625-55c0-4516-8ad2-d23fcd8ef345
mat(ConstGate.P1)

# ╔═╡ b2506c8a-618b-4c82-8f34-f3ee1dad3159
vizcircuit(control(3, 3, 1=>shift(π/4)))

# ╔═╡ f7b9023f-c019-49aa-ba7f-24c5ccf327ee
mat(control(3, 3, 1=>shift(π/4)))

# ╔═╡ eac02385-a859-4ac7-8fa3-c52798e673af
md"**Yao blocks are fast**"

# ╔═╡ 799fbfd4-f4a4-47ff-a7c1-9ccbbc4d5de4
let
	x = @bind run_benchmark_mat CheckBox()
	md" $x run the benchmark"
end

# ╔═╡ 8e517653-4657-4bbf-b606-956276d3b45f
if run_benchmark_mat
	@benchmark mat(control(16, 7, 3=>shift(π/4)))
end

# ╔═╡ 94c62812-6531-4b62-b9fa-6bbdcf37be56
if run_benchmark_mat
	@benchmark kron(
		kron(
			kron(Diagonal(ones(ComplexF64, 2^9)), Diagonal(mat(ConstGate.P1))),
				Diagonal(ones(ComplexF64, 2^3))
		), 
		kron(mat(ComplexF64, shift(π/4)), Diagonal(ones(ComplexF64, 2^2)))
	) +
	kron(
		kron(Diagonal(ones(ComplexF64, 2^9)), Diagonal(mat(ConstGate.P0))),
			Diagonal(ones(ComplexF64, 2^6))
	)
end

# ╔═╡ 97893414-d701-487e-9af1-09c06279c209
md"**Yao block is general purposed**"

# ╔═╡ 9be2d678-5557-471c-a78d-a3d38552d384
multi_control = control(10, (2, -6, 4), (9, 1)=>chain(put(2, 1=>X), rot(kron(X, X), 0.3)))

# ╔═╡ 1213dcc2-535b-4a1a-9582-fce6ef1ae556
vizcircuit(multi_control)

# ╔═╡ b3508525-ed3a-41f1-8ab8-6ebf2e055ac5
mat(multi_control)

# ╔═╡ 2d1d8978-90f4-4615-9d70-9c885154c269
md"""
`ChainBlock` is a composite block for matrix multiplication (reversed order), it is defined as:
```math
\texttt{chain(G₁, G₂, ..., Gₘ)} = G_m G_{m-1}\ldots G_1

```
"""

# ╔═╡ 1c2f0105-f482-488c-bafb-c1405c6ca3ee
qft3_first = chain(put(3, 1=>H), control(3, 2, 1=>shift(π/2)), control(3, 3, 1=>shift(π/4)))

# ╔═╡ 392a37a5-2ca2-4c58-88e9-a6eafcbdd622
vizcircuit(qft3_first)

# ╔═╡ 5256b92a-7a5b-4faf-9c0c-d1d8ace35b65
mat(qft3_first)

# ╔═╡ 142859a0-0dd1-4481-a250-033382ce71c4
mat(qft3)

# ╔═╡ e2020a46-03ac-4459-b7aa-8872ecb8a859
md"### Use matrix as a block"

# ╔═╡ 5fc34d66-bf11-4f18-80eb-ddc21e313830
random_gate = matblock(rand_unitary(2); tag="random gate")

# ╔═╡ 5a13c9ea-f20b-4b17-acf3-34502a795ee5
mat(put(3, 2=>random_gate))

# ╔═╡ 1c38da2b-411d-489e-96c1-ee43e66d60aa
vizcircuit(put(3, 2=>random_gate))

# ╔═╡ f1ca7769-5f45-49ab-8811-5ece0de21a58
md"### Blocks for d-level systems"

# ╔═╡ 247ce329-9960-49f6-812f-d537923606f9
block_3level = matblock(rand_unitary(3); nlevel=3, tag="3 level")

# ╔═╡ bf635286-f17a-43ad-af6d-f1dca4855ff8
mat(put(4, 2=>block_3level))

# ╔═╡ 12180746-4e54-4857-9fb3-60f60b8772c2
vizcircuit(put(4, 2=>block_3level))

# ╔═╡ 35bd48f3-a871-429a-bcce-fa7b1837c9cf
md"### Block properties"

# ╔═╡ 19247c86-4093-4e01-87f4-249b98938ef9
mat(X)

# ╔═╡ 58f78c85-b403-4b61-a228-0876785aa662
md"1. Hermitian: ``A^\dagger = A``"

# ╔═╡ 932f7b94-32ab-4f4a-9dcf-32e7acd52b1d
ishermitian(X)

# ╔═╡ cbc047b1-c2fa-4777-bb24-3361addd25e6
md"Only Hermitian operators are allowed in time evolution."

# ╔═╡ 16b49967-d21f-49f3-bdcc-ea0a598e2491
md"""
```math
\texttt{time\_evolve(H, t)} = e^{-iHt}
```
"""

# ╔═╡ 14c9e171-9a55-4967-a35c-dae6cebaa509
time_evolve(ConstGate.S, 0.4)

# ╔═╡ 3c186694-ffa1-490c-9009-74ed2916ce40
mat(ConstGate.S)  # this is not hermitian

# ╔═╡ 02b9e54f-30ff-4c21-a156-6dc18ea220e1
md"Hermitian operators are not changed by conjugate transpose"

# ╔═╡ de970766-3803-46d7-9d7c-8a4c5289cb15
qft3'

# ╔═╡ 29c7a849-fbed-4042-9e7b-2933553d0c01
vizcircuit(qft3')

# ╔═╡ 051ed1d9-25b2-46be-b484-2310d7e6af30
md"We can verify the correctness by checking the operator fidelity"

# ╔═╡ e8ffe842-eaa0-4a67-907f-21fe16cfe81a
operator_fidelity(qft3' * qft3, igate(3))

# ╔═╡ 4fbc0fdc-6eca-4f07-9853-52d9a5c29bb4
md"Operator fidelity is defined as
```math
F(A, B) = \frac{1}{D}|{\rm Tr}(A, B)|
```
where ``D`` is the Hiltert space dimension.
"

# ╔═╡ b6989319-f6aa-4d1c-b9cc-e4438a2eef5c
md"2. Reflexive: ``A^2 = \mathbb{1}``"

# ╔═╡ 8db5df39-1534-4752-b3af-b84527a4dd43
isreflexive(X)

# ╔═╡ 4b3ef73f-c8b4-43a1-8d98-2346a2d67f1a
md"""
```math
{\rm rot}(G, \theta) = \cos \frac{\theta}{2}I - i \sin \frac{\theta}{2} G
```
"""

# ╔═╡ a75b23bc-d3e1-4dba-bede-1f59b0a45af5
md"``\exp(iG\theta/2) = \uparrow``"

# ╔═╡ c3df22e6-fa28-48e4-bd6e-cab00099d44a
mat(rot(X, 0.4))

# ╔═╡ d4c491cd-de37-4c49-bbc7-d4d5db17374a
md"The Heisenberg Hamiltonian is hermitian"

# ╔═╡ e4b86d19-3d89-4410-99cb-d79feaa96e11
md"Only reflexive operators are allowed in rotation gate as the generator."

# ╔═╡ a04b4e93-d08a-4105-b27c-d5bee49296ff
rot(X + Y, 0.5)

# ╔═╡ ff0d1a92-5d70-41e7-a64a-65dc1e9b8e9a
X + Y

# ╔═╡ 05192a51-b578-43d3-a268-ab6bf4a5ffb9
md"3. Unitary: ``A A^\dagger = \mathbb{1}``"

# ╔═╡ 0b6d30ec-4ab1-4c96-9cd4-31cc30222483
isunitary(X)

# ╔═╡ 4c4925a2-2458-4423-858f-f13d819a207e
md"Unitary is required by the back propagation algorithm (will not be covered today)."

# ╔═╡ 04e1997e-ada3-4caa-8735-511422ab5c99
md"4. Commutative: ``[A, B] = \mathbb{0}``"

# ╔═╡ 04b12adc-b994-4cd1-8297-8b0fdc0a008b
iscommute(X, Y)

# ╔═╡ 7b527771-a8d2-4c0f-9987-3a033e4718e7
md"These properties does not require calling the `mat`, unless it falls back."

# ╔═╡ 40148613-9d2a-47b2-b696-f2999ec2cf14
iscommute(put(1000, 4=>X), put(1000, 6=>Y))

# ╔═╡ 09568a3a-d561-459c-adf4-cf49346f7294
iscommute(put(1000, 4=>X), put(1000, 4=>Y))

# ╔═╡ 98e69c16-63c3-4530-a480-59094f7fd0f6
ishermitian(EasyBuild.heisenberg(1000))

# ╔═╡ bda52396-f087-42f6-8c95-e3ba4b846f98
EasyBuild.heisenberg(3)

# ╔═╡ 2d43b615-463f-4717-856b-3bef30566622
md"## Run a circuit"

# ╔═╡ ef70a6af-56b0-4722-8b3e-039eb40227bc
md"### Run circuit on register"

# ╔═╡ 119c0675-939b-42be-8508-e077573373cb
reg = ghz_state(3)

# ╔═╡ 9ac72975-10bc-4f99-aced-5142c2931174
@with_terminal print_table(reg)

# ╔═╡ 33ddcca2-58af-4163-b6bd-139565d0b245
reg2 = apply(reg, qft3)

# ╔═╡ 117d2a17-296e-4252-a803-fad86b1496a3
md"Note: The inplace version `apply!` does not allocate and is faster, but Pluto does not track inplace variables very well."

# ╔═╡ c92bfbd4-2c72-4ce0-9bc7-aa671c05aa76
@with_terminal print_table(reg2)

# ╔═╡ 7ce66e78-d1d0-4b98-aa28-1ebce3cff27a
md"##### Run a subroutine"

# ╔═╡ 2ef0d491-f462-44fc-ae50-21901f770046
reg16 = rand_state(16)

# ╔═╡ d442292e-e6f2-49d3-a1f4-db2d72e27819
@time apply(reg16, put(16, 1:8=>EasyBuild.qft_circuit(8)))

# ╔═╡ d9e59f9a-4231-48b0-b28b-20fda13db4f2
@time apply(reg16, subroutine(16, EasyBuild.qft_circuit(8), 1:8))

# ╔═╡ ccf2576a-35a7-4ab3-95ac-c026fba6cc22
focus!(rand_state(10), (6,2,3)).state

# ╔═╡ f696c385-bcfd-416c-9078-fd74ff0d7a8b
md"A subroutine is similar to `put`, but runs faster when the sub-circuit size is large. This is because they use different algorithms."

# ╔═╡ 7c7e0ece-26a2-4f1b-a606-5d7d058b686b
md"### Run circuit on density matrix"

# ╔═╡ 7491bd3d-a212-43ba-8ad2-023468af6cb8
md"reduced density can be obtained by calling the `density_matrix` method."

# ╔═╡ d23b4fa0-82ba-4f8c-9aa8-5293b7421bfe
rho = density_matrix(reg, (2,3))

# ╔═╡ a1ba47de-699c-42a1-9d96-990163e21d1a
rho2 = apply(rho, put(2, 2=>X))

# ╔═╡ 8c12910d-8b97-4d00-956c-bf12253168ba
md"Density matrix supports unitary channel."

# ╔═╡ a8f7650e-2527-408f-9c44-ff0d67b52655
channel = unitary_channel([kron(X,Y), kron(Y,Z), kron(Z,X)], [0.4, 0.3, 0.3])

# ╔═╡ 317c4f26-861b-468e-9824-f1e0c3501b85
apply(rho, channel)

# ╔═╡ 8d7896a5-a452-48a0-8e37-1c58fff53fa8
md"""
Here, `kron` is a composite block that defined as
```math
\texttt{kron}(G_1, G_2, \ldots, G_n) = G_n \otimes G_{n-1} \otimes \ldots \otimes G_1
```
A `unitary_channel` does not have a matrix representation. When applying a unitary channel `unitary_channel(U, p)` to a density matrix, it effectively does
```math
\phi(\rho) = \sum_i p_i U_i ρ U_i^\dagger
```
"""

# ╔═╡ d861068a-5f1c-4800-9f02-f9d61137f30e
md"## Manipulate the block tree"

# ╔═╡ cb641f94-671d-4854-a45b-83c83c9dca07
md"### Recurse over the block tree"

# ╔═╡ 9683320e-07e4-4077-944f-cbbed4e13219
md"""Yao's block representation (QBIR) defines a tree, its siblings can be accessed with the `subblocks` function."""

# ╔═╡ 2858626d-c734-473e-b21d-66a818a87fb8
subblocks(qft3)

# ╔═╡ 21046c27-1f21-4fde-8376-8c0567be461c
subblocks(subblocks(qft3)[1])

# ╔═╡ 45b73eb4-b069-4c9e-906b-4e24702e2186
subblocks(subblocks(subblocks(qft3)[1])[1])

# ╔═╡ 254f7bc4-8c78-41ed-a7cd-cb24a1ad1397
subblocks(subblocks(subblocks(subblocks(qft3)[1])[1])[1])

# ╔═╡ cc3f7217-6f9a-43da-9124-e221d5195ba9
begin
	decompose(x::HGate)= Rz(0.5π)*Rx(0.5π)*Rz(0.5π)
	decompose(x::AbstractBlock)= chsubblocks(x, decompose.(subblocks(x)))
	vizcircuit(decompose(qft3); scale=0.65)
end

# ╔═╡ 9087fb14-8c2e-4d37-a6c7-12131536921a
md"Some methods are implemented by searching this tree in depth first order"

# ╔═╡ ad0bdf91-bf3d-4c82-9532-cd880cec3665
gatecount(qft3)

# ╔═╡ c2c0db82-2b70-41c1-8b04-a11fa3c964d0
parameters(qft3)

# ╔═╡ d15f73a7-75c2-4d60-8007-c6b492f9b7cd
dispatch(qft3, [1.0, 2.0, 3.0])

# ╔═╡ f12696ca-4739-4081-ab64-2ac1a8c83fe7
md"### Simplifying the block tree"

# ╔═╡ 0801fd0a-7216-44b8-a5fd-820db4b742e5
md"""
If you feel hard to handle two many types, you can compile away non-essential composite blocks like `repeat` and `kron` to basic types.
"""

# ╔═╡ fd2049d0-bf49-4f91-8198-fabf948e112a
Optimise.to_basictypes(kron(6, 1=>X, (2,3)=>SWAP))

# ╔═╡ 3d9d24ec-c7f1-45a1-833c-52942079f290
kron(6, 1=>X, (2,3)=>SWAP)

# ╔═╡ 790f16b5-ffab-45c9-b778-77f693d6bed5
Optimise.to_basictypes(repeat(10, X, (5,2)))

# ╔═╡ 50e4ab2e-8735-4176-8098-1955e1061e6a
md"However, these function does not propagate."

# ╔═╡ 8b209691-41e5-41bd-bb72-b06a0ff745ab
Optimise.to_basictypes(chain(repeat(10, X, (5,2))))

# ╔═╡ 15cf8e86-60b5-4a73-ac33-4317c2fc4b7b
md"One apply these functions as the simplification rules repeatedly through the tree"

# ╔═╡ b415cc69-5e54-42a3-9d2f-697fcc355173
Optimise.simplify(chain(repeat(10, X, (5,2))), rules=[Optimise.to_basictypes])

# ╔═╡ a1e1bf5a-6938-4942-8624-0b6021acf3eb
md"Apply this rule to the Heisenberg Hamiltonian"

# ╔═╡ 30891f23-c649-4fd7-b06a-68d18806c9d0
Optimise.simplify(EasyBuild.heisenberg(3), rules=[Optimise.to_basictypes])

# ╔═╡ 6d5ddcad-7e0f-4c5a-b740-0f11e60ff92d
md"You can also remove nested `ChainBlock` or `Add` block"

# ╔═╡ d24e0ff8-419f-4f14-ab64-aaf5ec02bbdf
Optimise.simplify(EasyBuild.heisenberg(3), rules=[Optimise.to_basictypes, Optimise.eliminate_nested])

# ╔═╡ 7ab63a69-1415-4772-adf1-d2ad95ea6323
md"## Play with Hamiltonians"

# ╔═╡ 9661f132-906e-4d48-8f11-00e91d0421bd
h30 = EasyBuild.heisenberg(30)

# ╔═╡ f371af51-c67e-4d19-90b9-abe3abb8214d
md"### Block indexing"

# ╔═╡ 6ef05e08-a91f-4a38-81da-09fc11826687
# j = DitStr{2}(rand(0:1, 30))  # programming way
# j = dit"0000000000_0000100000_0000000000;2"  # dit string
j = bit"0000000000_0000100000_0000000000"

# ╔═╡ 16e1e537-e00f-4d01-ac4d-ba088d0c86a1
h30[j, j]

# ╔═╡ 4f2835ec-1c06-4ba1-be9d-39866878457b
h30[:, j]

# ╔═╡ 9ad3aa3e-5aa3-4ec3-963d-ea1530e95859
h30[:, h30[:, j]]

# ╔═╡ d9f499ec-35e8-44c9-9124-20257ab3126e
let
	gadget = @bind num_propagate Slider(0:20; show_value=true)
	md"""num_propagate $gadget"""
end

# ╔═╡ c73b22d1-a4df-44cb-bbbc-a2b89c41bd06
let J = EntryTable([j], [1.0+0im])
	for i=1:num_propagate
		J = h30[:,J]
	end
	J
end

# ╔═╡ e0402159-f24a-48a4-af52-97ed545207e6
md"### The ground state"

# ╔═╡ 91c4e3a5-2ae2-4653-984f-1d6c802e04cb
md"### Time Evolution"

# ╔═╡ 6553204d-9e51-4cfa-8765-35d446889c9d
h15 = EasyBuild.heisenberg(15);

# ╔═╡ 2f38c792-aa13-41bd-a43e-f5bfdd9978f3
@time mat(h15)

# ╔═╡ cb934cb3-7703-4637-8c7a-dd5d027c6463
@time KrylovKit.eigsolve(mat(h15), statevec(rand_state(15)), 1, :SR)

# ╔═╡ 5101fa6a-f27b-405b-ba77-7468e558166c
@time let
	reg = rand_state(15)
	apply(reg, time_evolve(h15, 0.5))
end

# ╔═╡ 89cb4078-c993-4b26-81ce-2551f8f3c843
md"For larger Hamiltonian, we prefer storing the Hamiltonian matrix in the memory instead of computing the terms on the fly."

# ╔═╡ 2b38ee10-7110-4306-810e-969c7bfdc0c0
cache(h15)

# ╔═╡ d3ce1536-bcb5-4a76-b1f3-a7a7aa50dc9a
@time let
	reg = rand_state(15)
	apply(reg, time_evolve(cache(h15), 0.5))
end

# ╔═╡ 6d9263d5-0ef7-4ee9-a045-7324c4f2a5c9
md"## Eigen basis and measurement"

# ╔═╡ 4c6b0ade-7662-4c17-8d29-63d3e6f67371
YaoBlocks.eigenbasis(kron(4, 4=>X, 1=>X))

# ╔═╡ 89137cf0-be2a-417c-a7fb-9b87b01b8601
let
	reg = rand_state(4)
	res = measure!(reg)
	print_table(reg)
	res
end

# ╔═╡ b438d59a-b124-4d64-b006-53492512ee04
let
	reg = rand_state(4)
	res = measure!(kron(4, 4=>X, 1=>X), reg)
	print_table(apply(reg, kron(4, 4=>H, 1=>H)))  # print on the eigenbasis
	res  # the state projected to that eigenbasis
end

# ╔═╡ 70a83f9a-ee22-48da-9294-49e17992531f
md"## Extend the block system"

# ╔═╡ 04c0a186-0466-45ec-8185-82404cc83851
md"""
**General purposed**
* `YaoAPI.mat(T, block) -> AbstractMatrix`, matrix representation
* `YaoBlocks.unsafe_apply!(reg, block)`, allocation free method to apply a gate (optional)
* `YaoBlocks.print_block(io, block)`, better printing (optional)
**Composite blocks**
* `YaoAPI.occupied_locs(block) -> tuple`
* `YaoAPI.subblocks(block) -> tuple or vector`
* `YaoAPI.chsubblocks(block, subblocks)`

**Parameter management API**
* `setiparams!(block, parameters)`, setter of intrinsic parameters
* `getiparams(block) -> tuple`, getter of intrinsic parameters
"""

# ╔═╡ 2ba6c570-fc5e-49a1-b4cd-460b6ab625e0
md"""
An exmaple:
[https://github.com/QuantumBFS/Yao.jl/blob/master/src/EasyBuild/block_extension/FSimGate.jl](https://github.com/QuantumBFS/Yao.jl/blob/master/src/EasyBuild/block_extension/FSimGate.jl)
"""

# ╔═╡ 4eaa666d-a269-4513-bb28-0c55367e8921
md"## News
* Aqua test: [https://github.com/QuantumBFS/Yao.jl/issues/381](https://github.com/QuantumBFS/Yao.jl/issues/381)
* Next time: Bloqade [https://github.com/QuEraComputing/Bloqade.jl](https://github.com/QuEraComputing/Bloqade.jl)
* Announce of Bounty issue: [https://github.com/QuantumBFS/Yao.jl/issues/403](https://github.com/QuantumBFS/Yao.jl/issues/403)
"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
KrylovKit = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Yao = "5872b779-8223-5990-8dd0-5abbb0748c8c"
YaoPlots = "32cfe2d9-419e-45f2-8191-2267705d8dbc"

[compat]
BenchmarkTools = "~1.3.1"
KrylovKit = "~0.5.4"
PlutoUI = "~0.7.39"
Yao = "~0.8.1"
YaoPlots = "~0.7.4"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.0-rc1"
manifest_format = "2.0"
project_hash = "97ca8d6f76931008405ec0642f5c87492c4f126a"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "d0f59ebfe8d3ea2799fb3fb88742d69978e5843e"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.10"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "4c10eee4af024676200bc7752e536f858c6b8f93"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.1"

[[deps.BitBasis]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "2135f6117a102ef8ff951d01c3a826ec29c2b688"
uuid = "50ba71b6-fa0f-514d-ae9a-0916efc90dcf"
version = "0.8.0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CacheServers]]
deps = ["Distributed", "Test"]
git-tree-sha1 = "b584b04f236d3677b4334fab095796a128445bf8"
uuid = "a921213e-d44a-5460-ac04-5d720a99ba71"
version = "0.2.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9489214b993cd42d17f44c36e359bf6a7c919abf"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "1e315e3f4b0b7ce40feded39c73049692126cf53"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.3"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "0f4e115f6f34bbe43c19751c90a38b2f380637b9"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.3"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "9be8be1d8a6f44b96482c8af52238ea7987da3e3"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.45.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.Compose]]
deps = ["Base64", "Colors", "DataStructures", "Dates", "IterTools", "JSON", "LinearAlgebra", "Measures", "Printf", "Random", "Requires", "Statistics", "UUIDs"]
git-tree-sha1 = "9a2695195199f4f20b94898c8a8ac72609e165a4"
uuid = "a81c6b42-2e10-5240-aca2-a61377ecd94b"
version = "0.9.3"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Dierckx]]
deps = ["Dierckx_jll"]
git-tree-sha1 = "633c119fcfddf61fb4c75d77ce3ebab552a44723"
uuid = "39dd38d3-220a-591b-8e3c-4c3a8c710a94"
version = "0.5.2"

[[deps.Dierckx_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6596b96fe1caff3db36415eeb6e9d3b50bfe40ee"
uuid = "cd4c43a9-7502-52ba-aa6d-59fb2a88580b"
version = "0.1.0+0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.ExponentialUtilities]]
deps = ["ArrayInterfaceCore", "GPUArrays", "GenericSchur", "LinearAlgebra", "Printf", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "343c0b28b7513bbdd8ea91d8500fd1f357944f22"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.17.1"

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

[[deps.GPUArrays]]
deps = ["Adapt", "LLVM", "LinearAlgebra", "Printf", "Random", "Serialization", "Statistics"]
git-tree-sha1 = "c783e8883028bf26fb05ed4022c450ef44edd875"
uuid = "0c68f7d7-f131-5f86-a1c3-88cf8149b2d7"
version = "8.3.2"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "fb69b2a645fa69ba5f474af09221b9308b160ce6"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.3"

[[deps.GraphPlot]]
deps = ["ArnoldiMethod", "ColorTypes", "Colors", "Compose", "DelimitedFiles", "Graphs", "LinearAlgebra", "Random", "SparseArrays"]
git-tree-sha1 = "5cd479730a0cb01f880eff119e9803c13f214cab"
uuid = "a2cc645c-3eea-5389-862e-a155d0052231"
version = "0.5.2"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "4888af84657011a65afc7a564918d281612f983a"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.7.0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "b3364212fb5d870f724876ffcd34dd8ec6d98918"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.7"

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
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.KrylovKit]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "49b0c1dd5c292870577b8f58c51072bd558febb9"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.5.4"

[[deps.LLVM]]
deps = ["CEnum", "LLVMExtra_jll", "Libdl", "Printf", "Unicode"]
git-tree-sha1 = "e7e9184b0bf0158ac4e4aa9daf00041b5909bf1a"
uuid = "929cbde3-209d-540e-8aea-75f648917ca0"
version = "4.14.0"

[[deps.LLVMExtra_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg", "TOML"]
git-tree-sha1 = "771bfe376249626d3ca12bcd58ba243d3f961576"
uuid = "dad2f222-ce93-54a1-a47d-0025e8a3acab"
version = "0.0.16+0"

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

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "09e4b894ce6a976c354a69041a04748180d43637"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.15"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LuxurySparse]]
deps = ["InteractiveUtils", "LinearAlgebra", "Random", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "2ba5b1c20266ff288b0e3ecf718a77be5716ca90"
uuid = "d05aeea4-b7d4-55ac-b691-9e7fabb07ba2"
version = "0.6.13"

[[deps.MLStyle]]
git-tree-sha1 = "2041c1fd6833b3720d363c3ea8140bffaf86d9c4"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.12"

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
version = "2.28.0+0"

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
version = "2022.2.1"

[[deps.Multigraphs]]
deps = ["Graphs", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "055a7c49a626e17a8c99bcaaf472d0de60848929"
uuid = "7ebac608-6c66-46e6-9856-b5f43e107bac"
version = "0.3.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "1285416549ccfcdf0c50d4997a94331e88d68413"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8d1f54886b9037091edf146b517989fc4a09efec"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.39"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

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

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "383a578bdf6e6721f480e749d503ebc8405a0b22"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.6"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "2c11d7290036fe7aac9038ff312d3b3a2a5bf89e"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.4.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

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

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

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
deps = ["BitBasis", "LinearAlgebra", "LuxurySparse", "Reexport", "YaoAPI", "YaoArrayRegister", "YaoBlocks", "YaoSym"]
git-tree-sha1 = "d0f9e768783cc473c4c25eb753eceb418051a524"
uuid = "5872b779-8223-5990-8dd0-5abbb0748c8c"
version = "0.8.1"

[[deps.YaoAPI]]
git-tree-sha1 = "98cfa822c093a5727f0c9f6b94b9541e0c64a68c"
uuid = "0843a435-28de-4971-9e8b-a9641b2983a8"
version = "0.4.0"

[[deps.YaoArrayRegister]]
deps = ["Adapt", "BitBasis", "LegibleLambdas", "LinearAlgebra", "LuxurySparse", "MLStyle", "Random", "SparseArrays", "StaticArrays", "StatsBase", "TupleTools", "YaoAPI"]
git-tree-sha1 = "81ef9b5f217d960c5663091176631de5b6c3f7ec"
uuid = "e600142f-9330-5003-8abb-0ebd767abc51"
version = "0.9.1"

[[deps.YaoBlocks]]
deps = ["BitBasis", "CacheServers", "ChainRulesCore", "ExponentialUtilities", "InteractiveUtils", "LegibleLambdas", "LinearAlgebra", "LuxurySparse", "MLStyle", "Random", "SparseArrays", "StaticArrays", "StatsBase", "TupleTools", "YaoAPI", "YaoArrayRegister"]
git-tree-sha1 = "3dc8e5a7ebec4cc7f2854f8a66b5e77c5fc36fc1"
uuid = "418bc28f-b43b-5e0b-a6e7-61bbc1a2c1df"
version = "0.13.1"

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
deps = ["BitBasis", "Colors", "Compose", "GraphPlot", "Graphs", "Multigraphs", "Viznet", "Yao", "ZXCalculus"]
git-tree-sha1 = "49dd2db39a9c6e28d8d9c732c34619835c44c14c"
uuid = "32cfe2d9-419e-45f2-8191-2267705d8dbc"
version = "0.7.4"

[[deps.YaoSym]]
deps = ["BitBasis", "LinearAlgebra", "LuxurySparse", "Requires", "SparseArrays", "YaoArrayRegister", "YaoBlocks"]
git-tree-sha1 = "3518fa5611d252971771f4f4c015341f04367c89"
uuid = "3b27209a-d3d6-11e9-3c0f-41eb92b2cb9d"
version = "0.6.0"

[[deps.ZXCalculus]]
deps = ["Graphs", "LinearAlgebra", "MLStyle", "Multigraphs", "SparseArrays", "YaoHIR", "YaoLocations"]
git-tree-sha1 = "58e4f9a72618f2daf483f328fd82f0d10df8dc37"
uuid = "3525faa3-032d-4235-a8d4-8c2939a218dd"
version = "0.5.0"

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
# ╟─7eedb3d5-a3b6-40e5-8d2d-94623faf8d92
# ╟─1c835e36-299d-41c7-940f-6e8d9fea87f9
# ╟─84705dc3-dd2c-419b-983e-e6c879dd79da
# ╠═173125da-e761-11ec-25b9-fd328db9847e
# ╟─caec900e-13fd-4340-8c18-5c1d4ee91a29
# ╠═6b5100df-e813-4ed8-8f7c-c4fa05015ea6
# ╠═2d7486af-c5e7-465a-b8d9-d8f7b00d4db4
# ╠═c63152ea-9a68-4cff-bacf-978ef95ea9a8
# ╟─31619f1a-784b-48c4-9755-ff41d273f46f
# ╟─23a0f167-481e-4de4-b175-b96ee1923d7d
# ╠═63e8fea9-01c7-4f78-a493-acff4f412820
# ╟─0afe664a-2ce5-4867-b6c2-8169d492c4bc
# ╠═f2391ac6-8bc3-40c3-bb08-4604e4fd8dda
# ╠═221f1d4d-814a-4d51-9aa5-1a963acebf23
# ╟─f274ba5b-6a4a-40bf-af55-1a0bf10119a4
# ╠═d53ef0be-74c5-4f0e-820f-62ef4455f340
# ╟─3a739d22-f847-4459-bf1f-f7a58817a0ae
# ╟─729a8205-657d-4b1b-92a3-750ffe4ececb
# ╟─de0ea7f9-cd37-4e36-9bf3-dcc18658ad74
# ╠═0b4c480d-8c72-4437-b18b-7c87a3abdf32
# ╟─05063b40-e87f-49d8-ac15-66371b8e8cce
# ╠═ef0f8625-55c0-4516-8ad2-d23fcd8ef345
# ╠═b2506c8a-618b-4c82-8f34-f3ee1dad3159
# ╠═f7b9023f-c019-49aa-ba7f-24c5ccf327ee
# ╟─eac02385-a859-4ac7-8fa3-c52798e673af
# ╟─799fbfd4-f4a4-47ff-a7c1-9ccbbc4d5de4
# ╠═8e517653-4657-4bbf-b606-956276d3b45f
# ╠═94c62812-6531-4b62-b9fa-6bbdcf37be56
# ╟─97893414-d701-487e-9af1-09c06279c209
# ╠═9be2d678-5557-471c-a78d-a3d38552d384
# ╠═1213dcc2-535b-4a1a-9582-fce6ef1ae556
# ╠═b3508525-ed3a-41f1-8ab8-6ebf2e055ac5
# ╟─2d1d8978-90f4-4615-9d70-9c885154c269
# ╠═1c2f0105-f482-488c-bafb-c1405c6ca3ee
# ╠═392a37a5-2ca2-4c58-88e9-a6eafcbdd622
# ╠═5256b92a-7a5b-4faf-9c0c-d1d8ace35b65
# ╠═142859a0-0dd1-4481-a250-033382ce71c4
# ╟─e2020a46-03ac-4459-b7aa-8872ecb8a859
# ╠═5fc34d66-bf11-4f18-80eb-ddc21e313830
# ╠═5a13c9ea-f20b-4b17-acf3-34502a795ee5
# ╠═1c38da2b-411d-489e-96c1-ee43e66d60aa
# ╟─f1ca7769-5f45-49ab-8811-5ece0de21a58
# ╠═247ce329-9960-49f6-812f-d537923606f9
# ╠═bf635286-f17a-43ad-af6d-f1dca4855ff8
# ╠═12180746-4e54-4857-9fb3-60f60b8772c2
# ╟─35bd48f3-a871-429a-bcce-fa7b1837c9cf
# ╠═19247c86-4093-4e01-87f4-249b98938ef9
# ╟─58f78c85-b403-4b61-a228-0876785aa662
# ╠═932f7b94-32ab-4f4a-9dcf-32e7acd52b1d
# ╟─cbc047b1-c2fa-4777-bb24-3361addd25e6
# ╟─16b49967-d21f-49f3-bdcc-ea0a598e2491
# ╠═14c9e171-9a55-4967-a35c-dae6cebaa509
# ╠═3c186694-ffa1-490c-9009-74ed2916ce40
# ╟─02b9e54f-30ff-4c21-a156-6dc18ea220e1
# ╠═de970766-3803-46d7-9d7c-8a4c5289cb15
# ╠═29c7a849-fbed-4042-9e7b-2933553d0c01
# ╟─051ed1d9-25b2-46be-b484-2310d7e6af30
# ╠═e8ffe842-eaa0-4a67-907f-21fe16cfe81a
# ╟─4fbc0fdc-6eca-4f07-9853-52d9a5c29bb4
# ╟─b6989319-f6aa-4d1c-b9cc-e4438a2eef5c
# ╠═8db5df39-1534-4752-b3af-b84527a4dd43
# ╟─4b3ef73f-c8b4-43a1-8d98-2346a2d67f1a
# ╟─a75b23bc-d3e1-4dba-bede-1f59b0a45af5
# ╠═c3df22e6-fa28-48e4-bd6e-cab00099d44a
# ╟─d4c491cd-de37-4c49-bbc7-d4d5db17374a
# ╟─e4b86d19-3d89-4410-99cb-d79feaa96e11
# ╠═a04b4e93-d08a-4105-b27c-d5bee49296ff
# ╠═ff0d1a92-5d70-41e7-a64a-65dc1e9b8e9a
# ╟─05192a51-b578-43d3-a268-ab6bf4a5ffb9
# ╠═0b6d30ec-4ab1-4c96-9cd4-31cc30222483
# ╟─4c4925a2-2458-4423-858f-f13d819a207e
# ╟─04e1997e-ada3-4caa-8735-511422ab5c99
# ╠═04b12adc-b994-4cd1-8297-8b0fdc0a008b
# ╟─7b527771-a8d2-4c0f-9987-3a033e4718e7
# ╠═40148613-9d2a-47b2-b696-f2999ec2cf14
# ╠═09568a3a-d561-459c-adf4-cf49346f7294
# ╠═98e69c16-63c3-4530-a480-59094f7fd0f6
# ╠═bda52396-f087-42f6-8c95-e3ba4b846f98
# ╟─2d43b615-463f-4717-856b-3bef30566622
# ╟─ef70a6af-56b0-4722-8b3e-039eb40227bc
# ╠═119c0675-939b-42be-8508-e077573373cb
# ╠═9ac72975-10bc-4f99-aced-5142c2931174
# ╠═33ddcca2-58af-4163-b6bd-139565d0b245
# ╟─117d2a17-296e-4252-a803-fad86b1496a3
# ╠═c92bfbd4-2c72-4ce0-9bc7-aa671c05aa76
# ╟─7ce66e78-d1d0-4b98-aa28-1ebce3cff27a
# ╠═2ef0d491-f462-44fc-ae50-21901f770046
# ╠═d442292e-e6f2-49d3-a1f4-db2d72e27819
# ╠═d9e59f9a-4231-48b0-b28b-20fda13db4f2
# ╠═ccf2576a-35a7-4ab3-95ac-c026fba6cc22
# ╟─f696c385-bcfd-416c-9078-fd74ff0d7a8b
# ╟─7c7e0ece-26a2-4f1b-a606-5d7d058b686b
# ╟─7491bd3d-a212-43ba-8ad2-023468af6cb8
# ╠═d23b4fa0-82ba-4f8c-9aa8-5293b7421bfe
# ╠═a1ba47de-699c-42a1-9d96-990163e21d1a
# ╟─8c12910d-8b97-4d00-956c-bf12253168ba
# ╠═a8f7650e-2527-408f-9c44-ff0d67b52655
# ╠═317c4f26-861b-468e-9824-f1e0c3501b85
# ╟─8d7896a5-a452-48a0-8e37-1c58fff53fa8
# ╟─d861068a-5f1c-4800-9f02-f9d61137f30e
# ╟─cb641f94-671d-4854-a45b-83c83c9dca07
# ╟─9683320e-07e4-4077-944f-cbbed4e13219
# ╠═2858626d-c734-473e-b21d-66a818a87fb8
# ╠═21046c27-1f21-4fde-8376-8c0567be461c
# ╠═45b73eb4-b069-4c9e-906b-4e24702e2186
# ╠═254f7bc4-8c78-41ed-a7cd-cb24a1ad1397
# ╠═cc3f7217-6f9a-43da-9124-e221d5195ba9
# ╟─9087fb14-8c2e-4d37-a6c7-12131536921a
# ╠═ad0bdf91-bf3d-4c82-9532-cd880cec3665
# ╠═c2c0db82-2b70-41c1-8b04-a11fa3c964d0
# ╠═d15f73a7-75c2-4d60-8007-c6b492f9b7cd
# ╟─f12696ca-4739-4081-ab64-2ac1a8c83fe7
# ╟─0801fd0a-7216-44b8-a5fd-820db4b742e5
# ╠═fd2049d0-bf49-4f91-8198-fabf948e112a
# ╠═3d9d24ec-c7f1-45a1-833c-52942079f290
# ╠═790f16b5-ffab-45c9-b778-77f693d6bed5
# ╟─50e4ab2e-8735-4176-8098-1955e1061e6a
# ╠═8b209691-41e5-41bd-bb72-b06a0ff745ab
# ╟─15cf8e86-60b5-4a73-ac33-4317c2fc4b7b
# ╠═b415cc69-5e54-42a3-9d2f-697fcc355173
# ╟─a1e1bf5a-6938-4942-8624-0b6021acf3eb
# ╠═30891f23-c649-4fd7-b06a-68d18806c9d0
# ╟─6d5ddcad-7e0f-4c5a-b740-0f11e60ff92d
# ╠═d24e0ff8-419f-4f14-ab64-aaf5ec02bbdf
# ╟─7ab63a69-1415-4772-adf1-d2ad95ea6323
# ╠═9661f132-906e-4d48-8f11-00e91d0421bd
# ╟─f371af51-c67e-4d19-90b9-abe3abb8214d
# ╠═6ef05e08-a91f-4a38-81da-09fc11826687
# ╠═16e1e537-e00f-4d01-ac4d-ba088d0c86a1
# ╠═4f2835ec-1c06-4ba1-be9d-39866878457b
# ╠═9ad3aa3e-5aa3-4ec3-963d-ea1530e95859
# ╟─d9f499ec-35e8-44c9-9124-20257ab3126e
# ╠═c73b22d1-a4df-44cb-bbbc-a2b89c41bd06
# ╟─e0402159-f24a-48a4-af52-97ed545207e6
# ╠═2f38c792-aa13-41bd-a43e-f5bfdd9978f3
# ╠═cb934cb3-7703-4637-8c7a-dd5d027c6463
# ╟─91c4e3a5-2ae2-4653-984f-1d6c802e04cb
# ╠═6553204d-9e51-4cfa-8765-35d446889c9d
# ╠═5101fa6a-f27b-405b-ba77-7468e558166c
# ╟─89cb4078-c993-4b26-81ce-2551f8f3c843
# ╠═2b38ee10-7110-4306-810e-969c7bfdc0c0
# ╠═d3ce1536-bcb5-4a76-b1f3-a7a7aa50dc9a
# ╟─6d9263d5-0ef7-4ee9-a045-7324c4f2a5c9
# ╠═4c6b0ade-7662-4c17-8d29-63d3e6f67371
# ╠═89137cf0-be2a-417c-a7fb-9b87b01b8601
# ╠═b438d59a-b124-4d64-b006-53492512ee04
# ╟─70a83f9a-ee22-48da-9294-49e17992531f
# ╟─04c0a186-0466-45ec-8185-82404cc83851
# ╟─2ba6c570-fc5e-49a1-b4cd-460b6ab625e0
# ╟─4eaa666d-a269-4513-bb28-0c55367e8921
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
