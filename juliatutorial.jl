### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ ee916ff8-c4f8-4dfb-83c5-12d1ab95f111
using Pkg

# ╔═╡ 7d242d2a-d190-4a11-b218-60650ba70533
using PlutoUI

# ╔═╡ 52c27043-31c2-4e90-b6a5-d858aa7056d4
using AbstractTrees

# ╔═╡ d5d44e77-934f-4f0c-af1b-d89f0778142d
using Yao

# ╔═╡ 713939c6-4fe6-11ed-3e49-6bcc498b82f2
md"""
# Julia 语言基础
"""

# ╔═╡ 0919dfcc-b344-4e4c-abfa-9c3914e2850b
md"## 一些帮助函数"

# ╔═╡ 012b69d8-6304-4e91-9c0f-07fe3ad9980f
AbstractTrees.children(x::Type) = subtypes(x)

# ╔═╡ d4a6f68e-b7da-4ca1-b43c-c2da7929cd3d
function print_type_tree(T; maxdepth=5)
	io = IOBuffer()
	AbstractTrees.print_tree(io, T; maxdepth)
	Text(String(take!(io)))
end

# ╔═╡ 26b30265-558b-49e7-b9f5-0b8af30c1273
pkg_registries = Pkg.Operations.Context().registries;

# ╔═╡ 922071fb-dac2-436e-a343-d0d22bd3c864
function AbstractTrees.children(uuid::Base.UUID)
    dep = get(Pkg.dependencies(), uuid, nothing)
    values(dep.dependencies)
end

# ╔═╡ d75c0427-12fe-4b2d-9bd1-b08f477966a6
function AbstractTrees.printnode(io::IO, uuid::Base.UUID)
    dep = get(Pkg.dependencies(), uuid, nothing)
	link = collect(Pkg.Operations.find_urls(pkg_registries, uuid))
	if length(link) > 0
    	print(io, "<a href=\"$(link[1])\">$(dep.name)</a> (v$(dep.version))")
	else
		print(io, "$(dep.name)")
	end
end

# ╔═╡ a3b55379-b79b-433d-99c6-e617856de0f1


# ╔═╡ 9bb41efb-2817-4258-af2b-1fe515b6007a
macro mermaid_str(str)
	return HTML("""<script src="https://cdn.bootcss.com/mermaid/8.14.0/mermaid.min.js"></script>
<script>
  // how to do it correctly?
  mermaid.init({
    noteMargin: 10
  }, ".someClass");
</script>

<div class="mermaid someClass">
  $str
</div>
""")
end

# ╔═╡ a9a9f06e-4737-4619-b497-f488ea25fdf3
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

# ╔═╡ fb09bc52-7282-44c9-b4c1-eb0b02c287df
md"""
## 内容
* Julia 的软件包管理
* Julia 的语言特性
* Julia 与高性能计算
"""

# ╔═╡ 8e7f15fd-ae65-4559-972a-2c9720ac1547
md"# Julia 是什么样的语言?"

# ╔═╡ 73ce1dff-a3ff-431b-9acb-7af6c00b35f6
md"""
## 源代码开放的现代高性能语言

Julia 于2012年在 MIT 诞生， 其的源代码被托管在 [Github 仓库](https://github.com/JuliaLang/julia) 中， 其软件协议为可商用的 MIT 协议。 不仅 Julia 语言如此， 大多 Julia 的软件包系统也依托 Github 管理， 其协议也大多为开源。目的是为了解决两语言问题

* 速度 (🐟): C, C++, Fortran
* 开发效率 (🐾): Python, Matlab, Mathematica

What about **Python + C++**?平台移植差和可维护性变差。

![](https://user-images.githubusercontent.com/6257240/200309092-6a138366-ac52-47e5-a010-47711612632b.png)
"""

# ╔═╡ 8ea2593c-2f93-47c1-aa7d-918c848f8bfb
md"""
## 特别之处

Julia 有很多特别之处，在此列举一个其中最重要的一点

$(leftright(
	html"<div style='font-size:40px; text-align: center; padding-right=50px'>Mutliple<br>多重</div>",
	html"<div style='font-size:40px; text-align: center'> Dispatch<br>派发</div>",
))


"""

# ╔═╡ fa90e2ac-0040-453b-948b-4b1e6a1999e2
function Base.:*(x::Integer, s::String)
	string(x) * s
end

# ╔═╡ 5c4131eb-0f24-454d-b0ab-9addbd8cc1e9
function Base.:*(x::Integer, s::Integer)
	x * string(s)
end

# ╔═╡ e23b935b-eab0-4256-9983-84fab6ed6632
function print_dependency_tree(pkg; maxdepth=5)
	io = IOBuffer()
	AbstractTrees.print_tree(io, Pkg.project().dependencies[string(pkg)]; maxdepth)
	HTML("<p style='font-family: Consolas; line-height: 1.2em; max-height: 300px;'>" * replace(String(take!(io)), "\n"=>"<br>") * "</p>")
end

# ╔═╡ 44e76470-48ad-4869-b2bf-1838d7f4f560
5 * "x"

# ╔═╡ 8551813f-c3ee-4285-820a-af81bbe1c888
@which 5 * "x"

# ╔═╡ b3f72d4b-9f1f-46fd-8145-212f96c320f8
methods(*)[end-10:end]

# ╔═╡ 1cf7d604-04d9-4b52-a087-aa73c52093e5
@which 5 * 5

# ╔═╡ 7c4957ef-46f8-4896-84f8-1835595f5ba5
md"""
例： 
"""

# ╔═╡ 915a6f21-1d94-4aed-aaa3-3a58a34264d3
md"""## 看教程之前
以下内容不会在教程中涉及，但是看本教程的基础。

1. 你需要[配置 Julia 语言环境](https://github.com/CodingThrust/CodingClub/blob/main/1.julia-setup.md)。

2. 你需要配置 [Pluto notebook](https://github.com/fonsp/Pluto.jl) 以在本地打开该教程， 您可以通过[此链接]()下载本教程到本地。

3. 您最好对 Git 和 [GitHub](https://github.com/) 有基本的了解， 以便理解 Julia 的软件包管理系统。
"""

# ╔═╡ a72f4263-b034-4aa8-8611-d53166cbb718
md"""
# Julia 的软件生态
"""

# ╔═╡ 216d9db3-2d4a-47ef-89c6-70edfdd7bd53
mermaid"""
graph TD;
A["安装包命令 pkg> add Yao"] --> B["从 GitHub 更新 registry (如 General)"] --> C["解析依赖关系与版本并生成 Manifest.toml 文件"] --> D["从 GitHub 找到对应的软件仓库"]
D --> E["下载对应软件包的版本并安装"]
"""

# ╔═╡ d1b9aa30-ac64-4653-95b9-ab8695fbf34b
md"以量子计算软件包 Yao 为例， 它的依赖关系可以非常复杂。"

# ╔═╡ e61c0433-58b0-46bf-956d-41caecd70316
print_dependency_tree(Yao; maxdepth=2)

# ╔═╡ abab335a-011d-4d07-bbd1-6811fa231e01
md"### 多重派发"

# ╔═╡ c196dcdf-94dd-48e4-88fd-6c7769dc087a
md"基于多重派发的，动态编译的语言"

# ╔═╡ bdc2b727-6f19-4472-9090-ce998de3363c
md"### 类型树"

# ╔═╡ 373b8ec5-19f1-4948-ad1f-9f589d1394c9
print_type_tree(Number)

# ╔═╡ 7d3f8447-0dc9-41b5-bb96-9555d45c23f6
md"### 深度剖析"

# ╔═╡ a08d00ac-e57d-4ba2-b9ae-77ff2647f0c1
md"## 100 行实现自己的分子动力学模拟"

# ╔═╡ 30df0765-c59f-4833-9a1a-3c1c559357df
md"primitive type and composite type"

# ╔═╡ 3f0c4f16-b6c1-41e1-b4ec-7a1895261d53
sizeof(Int)

# ╔═╡ 8cef1ad0-f1ff-41f8-988e-297d2b3a223b


# ╔═╡ c53bff5e-02d7-4741-a891-897f149acad0
md"### 首先， 你需要一个类型"

# ╔═╡ 6a4e3859-b33e-4a0f-9b60-57ffbcf477da
struct Atom{D, T}
	location::NTuple{D, T}

	function Atom(location::NTuple{D, T}) where {D, T}
		return new{D, T}(location)
	end
	function Atom(location1::T, locations::T...) where {T<:Number}
		return Atom((location1, locations...))
	end
end

# ╔═╡ 380daffe-a1ce-48a4-b56b-c7318a78028b
md"""
知识点： Julia 中一个类型由类型民称和类型参数构成。我们拿数组类型 `Array{N,T}` 来举例。
"""

# ╔═╡ 696ca377-2c04-43c0-b993-e978873732ef
dump(Array{3, Float64})

# ╔═╡ aa2c6763-c4b7-4fe8-8baa-42a76a3807b1
md"重载类型的显示方式， "

# ╔═╡ bc06ee9c-658b-40fe-b1fd-309204baf45d
Base.show(io::IO, m::Atom) = print(io, "Atom at $(m.location)")

# ╔═╡ 1724c0ae-53b7-4c0f-8f95-12b8a76c1aaa
Atom((3.0,4.0))

# ╔═╡ 0cf01f1c-b572-4341-beaf-13dea1cd14c6
Atom(3.0,4.0)

# ╔═╡ fa213a07-a1a9-48c8-a7f1-e18ada33893b
md"""知识点：
在 Julia 中， `$` 符号代表取值。
这里 `"Molecure at $(m.location)"` 等价于 `"Molecure at " * string(m.location)`， 因为 Julia 里面用 `*` 代表 String 的连接。
"""

# ╔═╡ 42d28e14-9f85-4319-8dc1-f31780e8b4e9
location(m::Atom) = m.location

# ╔═╡ fb8b1be3-62f7-4875-baf5-135ae8f7de06
distance(m1::Atom, m2::Atom) = sum(location(m1) .- location(m2))

# ╔═╡ 2d8de79a-b835-40c1-8554-cd505e11be03
md"问题: 定义函数去访问数据会不会让程序更慢呀？"

# ╔═╡ ccd91e92-3a62-4afb-9933-44a2e1742769
md"""首先定义库伦相互作用势能函数
```math
E_c = \frac{q_1q_2}{|\vec r_1 - \vec r_2|}
```
"""

# ╔═╡ 989144bb-1822-415f-8b89-c4a83ec3ca84
md"""
范德瓦耳斯力

```math
E_v = \frac{-A}{|\vec r_1 - \vec r_2|^6}
```
"""

# ╔═╡ a177bd10-3941-49b3-bd80-9db1c4597fb2
md"电子之间的斥力

```math
E_e = \frac{A'}{|\vec r_1 - \vec r_2|^{12}}
```
"

# ╔═╡ 8846ce25-defc-47ad-8490-595f5090bd8a
coulomb(atom1, atom2) = charge(atom1) * charge(atom2) / distance(atom1, atom2)

# ╔═╡ 34ffecd6-202d-46af-862c-0bf34524aa63
md"""## 资源
### 交流
* Julia slack
* Julia discourse
* JuliaCN discourse

### 学习
* 安装/升级 Julia, 配置 IDE
### 这个 notebook
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AbstractTrees = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Yao = "5872b779-8223-5990-8dd0-5abbb0748c8c"

[compat]
AbstractTrees = "~0.4.3"
PlutoUI = "~0.7.48"
Yao = "~0.8.5"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "df9118afe0aa8250f4af6b3ec6f0b50aa37bc7fd"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "52b3b436f8f73133d7bc3a6c71ee7ed6ab2ab754"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.3"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "195c5505521008abea5aee4f96930717958eac6f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.4.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "c46fb7dd1d8ca1d213ba25848a5ec4e47a1a1b08"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.26"

[[deps.ArrayInterfaceGPUArrays]]
deps = ["Adapt", "ArrayInterfaceCore", "GPUArraysCore", "LinearAlgebra"]
git-tree-sha1 = "fc114f550b93d4c79632c2ada2924635aabfa5ed"
uuid = "6ba088a2-8465-4c0a-af30-387133b534db"
version = "0.2.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitBasis]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "f51ef0fdfa5d8643fb1c12df3899940fc8cf2bf4"
uuid = "50ba71b6-fa0f-514d-ae9a-0916efc90dcf"
version = "0.8.1"

[[deps.CacheServers]]
deps = ["Distributed", "Test"]
git-tree-sha1 = "b584b04f236d3677b4334fab095796a128445bf8"
uuid = "a921213e-d44a-5460-ac04-5d720a99ba71"
version = "0.2.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e7ff6cadf743c098e08fca25c91103ee4303c9bb"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.6"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "3ca828fe1b75fa84b021a7860bd039eaea84d2f2"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.3.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.DataAPI]]
git-tree-sha1 = "46d2680e618f8abd007bce0c3026cb0c4a8f2032"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.12.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "c36550cb29cbe373e95b3f40486b9a4148f89ffd"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.2"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.ExponentialUtilities]]
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceGPUArrays", "GPUArraysCore", "GenericSchur", "LinearAlgebra", "Printf", "SparseArrays", "libblastrampoline_jll"]
git-tree-sha1 = "9837d3f3a904c7a7ab9337759c0093d3abea1d81"
uuid = "d4d017d3-3776-5f7e-afef-a10c40355c18"
version = "1.22.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "6872f5ec8fd1a38880f027a26739d42dcda6691f"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.2"

[[deps.GenericSchur]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "fb69b2a645fa69ba5f474af09221b9308b160ce6"
uuid = "c145ed77-6b09-5dd9-b285-bf645a82121e"
version = "0.5.3"

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

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

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
version = "7.84.0+0"

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
git-tree-sha1 = "94d9c52ca447e23eac0c0f074effbcd38830deb5"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.18"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LuxurySparse]]
deps = ["InteractiveUtils", "LinearAlgebra", "Random", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "660da52355791ea967982f86fd15aa8b4c9eae6d"
uuid = "d05aeea4-b7d4-55ac-b691-9e7fabb07ba2"
version = "0.7.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MLStyle]]
git-tree-sha1 = "43f9be9c281179fe44205e2dc19f22e71e022d41"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.15"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

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
git-tree-sha1 = "6c01a9b494f6d2a9fc180a08b182fcb06f0958a0"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.4.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "efc140104e6d0ae3e7e30d56c98c4a927154d684"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.48"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

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
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "f86b3a049e5d05227b10e15dbb315c5b90f14988"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.9"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

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
version = "1.10.1"

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

[[deps.URIs]]
git-tree-sha1 = "e59ecc5a41b000fa94423a578d29290c7266fc10"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Yao]]
deps = ["BitBasis", "LinearAlgebra", "LuxurySparse", "Reexport", "YaoAPI", "YaoArrayRegister", "YaoBlocks", "YaoSym"]
git-tree-sha1 = "58573a875eb3705c752de1ac3e4e228e7cfbc781"
uuid = "5872b779-8223-5990-8dd0-5abbb0748c8c"
version = "0.8.5"

[[deps.YaoAPI]]
git-tree-sha1 = "4732ed765411aef7983123961d34cd9e9729da4f"
uuid = "0843a435-28de-4971-9e8b-a9641b2983a8"
version = "0.4.3"

[[deps.YaoArrayRegister]]
deps = ["Adapt", "BitBasis", "DocStringExtensions", "LegibleLambdas", "LinearAlgebra", "LuxurySparse", "MLStyle", "Random", "SparseArrays", "StaticArrays", "StatsBase", "TupleTools", "YaoAPI"]
git-tree-sha1 = "ef1054c7d6dd71c184c068c04ce862f86f9a468b"
uuid = "e600142f-9330-5003-8abb-0ebd767abc51"
version = "0.9.3"

[[deps.YaoBlocks]]
deps = ["BitBasis", "CacheServers", "ChainRulesCore", "DocStringExtensions", "ExponentialUtilities", "InteractiveUtils", "LegibleLambdas", "LinearAlgebra", "LuxurySparse", "MLStyle", "Random", "SparseArrays", "StaticArrays", "StatsBase", "TupleTools", "YaoAPI", "YaoArrayRegister"]
git-tree-sha1 = "6d991dc024d604c2cdb6746ea71d8781c10b1a03"
uuid = "418bc28f-b43b-5e0b-a6e7-61bbc1a2c1df"
version = "0.13.5"

[[deps.YaoSym]]
deps = ["BitBasis", "LinearAlgebra", "LuxurySparse", "Requires", "SparseArrays", "YaoArrayRegister", "YaoBlocks"]
git-tree-sha1 = "118e2c434e810dd52a3564a1b99f7fd3a2bbb63e"
uuid = "3b27209a-d3d6-11e9-3c0f-41eb92b2cb9d"
version = "0.6.2"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╟─713939c6-4fe6-11ed-3e49-6bcc498b82f2
# ╟─0919dfcc-b344-4e4c-abfa-9c3914e2850b
# ╠═7d242d2a-d190-4a11-b218-60650ba70533
# ╠═52c27043-31c2-4e90-b6a5-d858aa7056d4
# ╠═012b69d8-6304-4e91-9c0f-07fe3ad9980f
# ╠═d4a6f68e-b7da-4ca1-b43c-c2da7929cd3d
# ╠═ee916ff8-c4f8-4dfb-83c5-12d1ab95f111
# ╠═26b30265-558b-49e7-b9f5-0b8af30c1273
# ╠═922071fb-dac2-436e-a343-d0d22bd3c864
# ╠═d75c0427-12fe-4b2d-9bd1-b08f477966a6
# ╠═e23b935b-eab0-4256-9983-84fab6ed6632
# ╠═a3b55379-b79b-433d-99c6-e617856de0f1
# ╠═9bb41efb-2817-4258-af2b-1fe515b6007a
# ╠═a9a9f06e-4737-4619-b497-f488ea25fdf3
# ╟─fb09bc52-7282-44c9-b4c1-eb0b02c287df
# ╟─8e7f15fd-ae65-4559-972a-2c9720ac1547
# ╟─73ce1dff-a3ff-431b-9acb-7af6c00b35f6
# ╟─8ea2593c-2f93-47c1-aa7d-918c848f8bfb
# ╠═fa90e2ac-0040-453b-948b-4b1e6a1999e2
# ╠═44e76470-48ad-4869-b2bf-1838d7f4f560
# ╠═8551813f-c3ee-4285-820a-af81bbe1c888
# ╠═b3f72d4b-9f1f-46fd-8145-212f96c320f8
# ╠═1cf7d604-04d9-4b52-a087-aa73c52093e5
# ╠═5c4131eb-0f24-454d-b0ab-9addbd8cc1e9
# ╠═7c4957ef-46f8-4896-84f8-1835595f5ba5
# ╟─915a6f21-1d94-4aed-aaa3-3a58a34264d3
# ╟─a72f4263-b034-4aa8-8611-d53166cbb718
# ╟─216d9db3-2d4a-47ef-89c6-70edfdd7bd53
# ╟─d1b9aa30-ac64-4653-95b9-ab8695fbf34b
# ╠═d5d44e77-934f-4f0c-af1b-d89f0778142d
# ╠═e61c0433-58b0-46bf-956d-41caecd70316
# ╟─abab335a-011d-4d07-bbd1-6811fa231e01
# ╟─c196dcdf-94dd-48e4-88fd-6c7769dc087a
# ╟─bdc2b727-6f19-4472-9090-ce998de3363c
# ╠═373b8ec5-19f1-4948-ad1f-9f589d1394c9
# ╟─7d3f8447-0dc9-41b5-bb96-9555d45c23f6
# ╟─a08d00ac-e57d-4ba2-b9ae-77ff2647f0c1
# ╟─30df0765-c59f-4833-9a1a-3c1c559357df
# ╠═3f0c4f16-b6c1-41e1-b4ec-7a1895261d53
# ╠═8cef1ad0-f1ff-41f8-988e-297d2b3a223b
# ╟─c53bff5e-02d7-4741-a891-897f149acad0
# ╠═6a4e3859-b33e-4a0f-9b60-57ffbcf477da
# ╟─380daffe-a1ce-48a4-b56b-c7318a78028b
# ╠═696ca377-2c04-43c0-b993-e978873732ef
# ╠═aa2c6763-c4b7-4fe8-8baa-42a76a3807b1
# ╠═bc06ee9c-658b-40fe-b1fd-309204baf45d
# ╠═1724c0ae-53b7-4c0f-8f95-12b8a76c1aaa
# ╠═0cf01f1c-b572-4341-beaf-13dea1cd14c6
# ╟─fa213a07-a1a9-48c8-a7f1-e18ada33893b
# ╠═42d28e14-9f85-4319-8dc1-f31780e8b4e9
# ╠═fb8b1be3-62f7-4875-baf5-135ae8f7de06
# ╟─2d8de79a-b835-40c1-8554-cd505e11be03
# ╟─ccd91e92-3a62-4afb-9933-44a2e1742769
# ╟─989144bb-1822-415f-8b89-c4a83ec3ca84
# ╟─a177bd10-3941-49b3-bd80-9db1c4597fb2
# ╠═8846ce25-defc-47ad-8490-595f5090bd8a
# ╠═34ffecd6-202d-46af-862c-0bf34524aa63
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
