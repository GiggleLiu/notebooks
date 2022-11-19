### A Pluto.jl notebook ###
# v0.19.15

using Markdown
using InteractiveUtils

# ╔═╡ ee916ff8-c4f8-4dfb-83c5-12d1ab95f111
using Pkg

# ╔═╡ 7d242d2a-d190-4a11-b218-60650ba70533
using PlutoUI

# ╔═╡ 52c27043-31c2-4e90-b6a5-d858aa7056d4
using AbstractTrees

# ╔═╡ cf0eb0cd-bcb7-4f7c-b462-bef13d3c2a97
using Libdl

# ╔═╡ 9954c036-d4d3-42c9-acbf-22623f84f254
using PyCall

# ╔═╡ c73baba2-9ec7-461e-b4e7-fd162606e134
using BenchmarkTools

# ╔═╡ d5d44e77-934f-4f0c-af1b-d89f0778142d
using Yao

# ╔═╡ 713939c6-4fe6-11ed-3e49-6bcc498b82f2
md"""
# 给 Julia 开发者的入门教程
"""

# ╔═╡ 0919dfcc-b344-4e4c-abfa-9c3914e2850b
md"## 一些帮助函数"

# ╔═╡ 156a1a62-e131-403f-b2a2-80f49e6a9b33
html"<button onclick=present()>Present</button>"

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

# ╔═╡ e23b935b-eab0-4256-9983-84fab6ed6632
function print_dependency_tree(pkg; maxdepth=5)
	io = IOBuffer()
	AbstractTrees.print_tree(io, Pkg.project().dependencies[string(pkg)]; maxdepth)
	HTML("<p style='font-family: Consolas; line-height: 1.2em; max-height: 300px;'>" * replace(String(take!(io)), "\n"=>"<br>") * "</p>")
end

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

# ╔═╡ bb346eb2-e070-4522-a991-1bfd0c2b05dc
function livecoding(src)
	HTML("""
<link rel="stylesheet" type="text/css" href="https://github.com/asciinema/asciinema-player/releases/download/v3.0.1/asciinema-player.css" />
<div id="demo"></div>
<script src="https://github.com/asciinema/asciinema-player/releases/download/v3.0.1/asciinema-player.min.js"></script>
<script>
AsciinemaPlayer.create('$src', document.getElementById('demo'));
</script>
""")
end

# ╔═╡ 4e054cdc-8004-4e2a-9eb4-e30b52992386
livecoding("https://raw.githubusercontent.com/GiggleLiu/YaoTutorial/master/clips/yao.cast")

# ╔═╡ fb09bc52-7282-44c9-b4c1-eb0b02c287df
md"""
## 内容
* Julia 语言简介
* Julia 的语言特性
* Julia 语言开发者
"""

# ╔═╡ 1ab95944-524b-43d8-a95e-da345634f4c1
md"""
[配置开发环境 - 中文版](https://discourse.juliacn.com/t/topic/6806)
"""

# ╔═╡ 8e7f15fd-ae65-4559-972a-2c9720ac1547
md"# Julia 是什么样的语言?"

# ╔═╡ 73ce1dff-a3ff-431b-9acb-7af6c00b35f6
md"""
## 源代码开放的现代高性能语言

Julia 于2012年在 MIT 诞生。Julia 语言的解释/编译器的源代码是开放的，被托管在 [Github 仓库](https://github.com/JuliaLang/julia) 中，其软件协议为可商用的 MIT 协议。 不仅 Julia 语言如此， 大多 Julia 的软件包系统也依托 Github 管理， 其协议也大多为开源。

Julia 语言被设计出来的目的是为了兼顾代码执行速度与开发效率.

* 执行速度: C, C++, Fortran
* 开发效率: Python, $(html"Matlab")
"""

# ╔═╡ ff0a8030-9a18-4d27-9a87-bed9aed0d2a8
md"# 编译语言快的秘诀"

# ╔═╡ fe174dbe-5c4b-4445-b485-5c21cc1e8917
md"静态类型程序的执行很快，因为类型信息被提前知道就可以被高效的编译。"

# ╔═╡ 000b93e6-8a1d-4c67-b5da-5013c6421e2c
mermaid"""
flowchart LR;
A("一段静态类型程序") --> | 编译/很慢 | B("二进制文件") --> | 执行/快 | C(结果)
"""

# ╔═╡ e4c3c93b-f2a7-4e0d-acb2-2a2d40b90385
const Clib = tempname()

# ╔═╡ 33a43668-4484-47d2-a7a6-09d930232252
let
	# prepare the source code
	source_name = "$Clib.c"
	open(source_name, "w") do f
		write(f, """
#include <stddef.h>
int c_factorial(size_t n) {
	int s = 1;
	for (size_t i=1; i<=n; i++) {
		s *= i;
	}
	return s;
}
""")
	end
	# compile to a shared library by piping C_code to gcc;
	# (only works if you have gcc installed)
	run(`gcc $source_name -fPIC -O3 -msse3 -shared -o $(Clib * "." * Libdl.dlext)`)
end

# ╔═╡ 2a22f131-6a99-4744-8914-19c8776700e7
c_factorial(x) = @ccall Clib.c_factorial(x::Csize_t)::Int

# ╔═╡ ab045ed0-7cbb-4565-bd7f-239dd94ce99e
md"# 解释语言方便的秘诀"

# ╔═╡ f3695873-435d-44cb-b9fb-af34dc38bdaa
md"动态类型的语言它不需要被编译"

# ╔═╡ ef736f15-6180-46ed-ac52-d57ac17429e8
mermaid"""
flowchart LR;
A("一段静态类型程序") --> | 解释执行/慢 | C(结果)
"""

# ╔═╡ 0f526702-f8e6-492d-bd14-e81874e6fefe
py"""
def factorial(n):
	x = 1
	for i in range(1, n+1):
	    x = x * i
	return x
"""

# ╔═╡ 46cf6881-650e-4ba1-a0dc-bcda67fb367b
py"factorial"(1000)

# ╔═╡ 922a2063-f516-46a5-95a9-9e0adca018aa
# `typemax` 可以获取类型的最大值
typemax(Int)

# ╔═╡ 105852eb-8f34-4d52-8ec3-68dff6997efb
md"🤔"

# ╔═╡ e6fd7a35-e45e-4cc7-ae24-7c2f8fd7c73d
md"但由于数据没有固定的类型，解释执行的语言必须用一个`Box(type, *data)`来表示一个数据。"

# ╔═╡ f7e5304d-7573-4e8c-b516-4c16a7432067
md"""## 双语言 **Python & C++** 的问题?
### 可维护性变差
* 配置 setup 文件更加复杂, 平台移植性变差，
* 培养新人成本过高,

### 非常适合张量运算, 但很多程序抽象发生在底层
* 蒙特卡洛 (Monte Carlo) 和模拟退火 (Simulated Annealing) 方法, 频繁多变的随机数生成和采样
* 范型张量网络 (Generic Tensor Network), 张量中的标量类型的基本元素非实数, 而是 tropical number 或者有限域 (Finite Field Algebra)
* 组合优化中的分支界定法 (branching)

![](https://user-images.githubusercontent.com/6257240/200309092-6a138366-ac52-47e5-a010-47711612632b.png)
"""

# ╔═╡ d04b2eca-9662-4518-8bb6-8b1bf07e8984
mermaid"""
flowchart LR;
A("一个 Julia 函数") --> B{有函数实例?}
B -- 否 --> N[推导数据类型<br>并编译/不快] --> C("内存中的二进制码")
C --> |执行/快| Z("结果")
B -- 是 --> C
"""

# ╔═╡ 3e3a2f23-8098-4d06-b4d1-157c97e4c094
md"函数实例 (method instance)： 内存中，一个针对特定输入类型的函数被编译后的二进制码。"

# ╔═╡ be4da897-df85-4276-bde1-7c1824cae796
md"""
### Julia 函数被编译的过程
0. 拿到一段 Julia 程序
1. 在 Julia 的中间表示 (Intermediate Representation) 上推导数据类型
2. 将带类型的程序编译到 LLVM 中间表示上
3. LLVM 在这个中间表示的基础上生成高性能的汇编/二进制码
"""

# ╔═╡ 04b5f8fc-32c1-430c-8bec-3e1a06bdda24


# ╔═╡ 13bcf3d6-2418-46e1-acde-050914064741
function jlfactorial(n)
	x = 1
	for i in 1:n
    	x = x * i
	end
	return x
end

# ╔═╡ 4253af25-41bd-47b6-a11e-c2902c677963
jlfactorial(1000)

# ╔═╡ db779958-e7d5-4164-87a7-219257ae45f0
@code_typed jlfactorial(1000)

# ╔═╡ 7b8e9026-6dc1-4d28-a2a7-912399a4fd51
with_terminal() do
	@code_native jlfactorial(1000)
end

# ╔═╡ 01972597-9d31-4972-a15d-51832f0f5910
@benchmark c_factorial(1000)

# ╔═╡ ec33aba5-28c9-4be9-9804-361f65de1f7a
@benchmark jlfactorial(1000)

# ╔═╡ 79e3c220-c281-4ab0-988a-39e1b0a39d64
@benchmark $(py"factorial")(1000)

# ╔═╡ 917e187d-5eda-49d6-a72a-0ed3f60d82d6
md"[learn more about calling C code](https://docs.julialang.org/en/v1/manual/calling-c-and-fortran-code/)"

# ╔═╡ 2613110d-1ab8-413a-85ce-a2985ee420da
md"## Julia 语言被编译到了 LLVM"

# ╔═╡ fb53a9ed-df58-410a-8275-e15718514950
md"""
LLVM 是很多语言的后端 Julia, Rust, Swift, Kotlin et al.
"""

# ╔═╡ e5b59cc9-0d14-4d8a-bb25-738540e7ebf9
with_terminal() do 
	@code_llvm jlfactorial(10)
end

# ╔═╡ ff27c9fc-0e55-47dc-b189-534f7a48fd3f
md"![LLVM](https://upload.wikimedia.org/wikipedia/en/d/dd/LLVM_logo.png)"

# ╔═╡ 8ea2593c-2f93-47c1-aa7d-918c848f8bfb
md"""
## 特别之处

Julia 有很多特别之处，在此列举一个其中最重要的一点

$(leftright(
	html"<div style='font-size:40px; text-align: center; padding-right=50px'>Mutliple<br>多重</div>",
	html"<div style='font-size:40px; text-align: center'> Dispatch<br>派发</div>",
))


"""

# ╔═╡ 69fed6cc-030b-4066-a023-0bbf1637fbbc
begin
	function mul(x::Real, s::Real)
		string(x) * " × " * string(s)
	end
	function f(x::Real, s::Real)
		mul(x, s)
	end
	function f(x::Float64, s::Real)
		mul(round(x; ndigits=3), s)
	end
	function f(x::Real, s::Float64)
		mul(x, round(s; ndigits=3))
	end
	function f(x::Real, s::Union{Float64, ComplexF64})
		mul(x, round(s; ndigits=3))
	end
end

# ╔═╡ 46cd1ee1-e269-46a7-93d3-72597b53a9a9
Float64 <: Real

# ╔═╡ 0e6cf095-5249-4660-8d92-1347b143f795
Float64 <: Union{Float64, ComplexF64}

# ╔═╡ 3290dfe4-0ec1-490a-92f9-a2b43d0ae344
Union{Float64, ComplexF64} <: Real

# ╔═╡ f3a267f5-f4cd-4056-a16a-538c2ae4d756
Union{Float64, ComplexF64} <: Number

# ╔═╡ c04c4d58-0469-45cc-a217-444a2b607245
print_type_tree(Real)

# ╔═╡ b3f72d4b-9f1f-46fd-8145-212f96c320f8
methods(f)

# ╔═╡ fed9cd8e-2eed-4195-8bab-f6c1b0d3785d
f(5, 6)

# ╔═╡ 1cf7d604-04d9-4b52-a087-aa73c52093e5
@which f(5, 5)

# ╔═╡ 5e276fd0-887e-4de2-b502-359be36e6fb6
md"最具体的获胜"

# ╔═╡ 4b12f0d9-e4e1-4214-9127-40612f38d7a3
@which f(5, Float64(π))

# ╔═╡ 9b00810e-8dc8-4602-a185-28e60c027b99
md"有时候，难论输赢"

# ╔═╡ 8c683b66-1fb2-49ad-9caf-cb891520f5c6
f(Float64(5), Float64(5))

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

# ╔═╡ daa6c126-c5ea-4f9c-b04d-39da54e3fc4c
html"""
<link rel="stylesheet" type="text/css" href="https://github.com/asciinema/asciinema-player/releases/download/v3.0.1/asciinema-player.css" />
<div id="demo"></div>
<script src="https://github.com/asciinema/asciinema-player/releases/download/v3.0.1/asciinema-player.min.js"></script>
<script>
AsciinemaPlayer.create('https://raw.githubusercontent.com/GiggleLiu/YaoTutorial/master/clips/yao.cast', document.getElementById('demo'));
</script>
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AbstractTrees = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
BenchmarkTools = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
Libdl = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
Yao = "5872b779-8223-5990-8dd0-5abbb0748c8c"

[compat]
AbstractTrees = "~0.4.3"
BenchmarkTools = "~1.3.2"
PlutoUI = "~0.7.48"
PyCall = "~1.94.1"
Yao = "~0.8.5"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "8a3ab74d63acc8d27d4daa56dcc44f532108c58f"

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

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "d9a9701b899b30332bbcb3e1679c41cce81fb0e8"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.2"

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

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "6e47d11ea2776bc5627421d59cdcc1296c058071"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.7.0"

[[deps.DataAPI]]
git-tree-sha1 = "e08915633fcb3ea83bf9d6126292e5bc5c739922"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.13.0"

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
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "cceb0257b662528ecdf0b4b4302eb00e767b38e7"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.0"

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

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "53b8b07b721b77144a0fbbbc2675222ebf40a02d"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.94.1"

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

[[deps.SnoopPrecompile]]
git-tree-sha1 = "f604441450a3c0569830946e5b33b78c928e1a85"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

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

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

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
# ╠═156a1a62-e131-403f-b2a2-80f49e6a9b33
# ╠═7d242d2a-d190-4a11-b218-60650ba70533
# ╠═52c27043-31c2-4e90-b6a5-d858aa7056d4
# ╠═012b69d8-6304-4e91-9c0f-07fe3ad9980f
# ╠═d4a6f68e-b7da-4ca1-b43c-c2da7929cd3d
# ╠═ee916ff8-c4f8-4dfb-83c5-12d1ab95f111
# ╠═26b30265-558b-49e7-b9f5-0b8af30c1273
# ╠═922071fb-dac2-436e-a343-d0d22bd3c864
# ╠═d75c0427-12fe-4b2d-9bd1-b08f477966a6
# ╠═e23b935b-eab0-4256-9983-84fab6ed6632
# ╠═9bb41efb-2817-4258-af2b-1fe515b6007a
# ╠═a9a9f06e-4737-4619-b497-f488ea25fdf3
# ╠═bb346eb2-e070-4522-a991-1bfd0c2b05dc
# ╠═4e054cdc-8004-4e2a-9eb4-e30b52992386
# ╟─fb09bc52-7282-44c9-b4c1-eb0b02c287df
# ╟─1ab95944-524b-43d8-a95e-da345634f4c1
# ╟─8e7f15fd-ae65-4559-972a-2c9720ac1547
# ╟─73ce1dff-a3ff-431b-9acb-7af6c00b35f6
# ╟─ff0a8030-9a18-4d27-9a87-bed9aed0d2a8
# ╟─fe174dbe-5c4b-4445-b485-5c21cc1e8917
# ╟─000b93e6-8a1d-4c67-b5da-5013c6421e2c
# ╠═e4c3c93b-f2a7-4e0d-acb2-2a2d40b90385
# ╠═cf0eb0cd-bcb7-4f7c-b462-bef13d3c2a97
# ╠═33a43668-4484-47d2-a7a6-09d930232252
# ╠═2a22f131-6a99-4744-8914-19c8776700e7
# ╟─ab045ed0-7cbb-4565-bd7f-239dd94ce99e
# ╟─f3695873-435d-44cb-b9fb-af34dc38bdaa
# ╟─ef736f15-6180-46ed-ac52-d57ac17429e8
# ╠═9954c036-d4d3-42c9-acbf-22623f84f254
# ╠═0f526702-f8e6-492d-bd14-e81874e6fefe
# ╠═46cf6881-650e-4ba1-a0dc-bcda67fb367b
# ╠═922a2063-f516-46a5-95a9-9e0adca018aa
# ╟─105852eb-8f34-4d52-8ec3-68dff6997efb
# ╟─e6fd7a35-e45e-4cc7-ae24-7c2f8fd7c73d
# ╟─f7e5304d-7573-4e8c-b516-4c16a7432067
# ╟─d04b2eca-9662-4518-8bb6-8b1bf07e8984
# ╟─3e3a2f23-8098-4d06-b4d1-157c97e4c094
# ╟─be4da897-df85-4276-bde1-7c1824cae796
# ╠═04b5f8fc-32c1-430c-8bec-3e1a06bdda24
# ╠═13bcf3d6-2418-46e1-acde-050914064741
# ╠═4253af25-41bd-47b6-a11e-c2902c677963
# ╠═db779958-e7d5-4164-87a7-219257ae45f0
# ╠═7b8e9026-6dc1-4d28-a2a7-912399a4fd51
# ╠═c73baba2-9ec7-461e-b4e7-fd162606e134
# ╠═01972597-9d31-4972-a15d-51832f0f5910
# ╠═ec33aba5-28c9-4be9-9804-361f65de1f7a
# ╠═79e3c220-c281-4ab0-988a-39e1b0a39d64
# ╟─917e187d-5eda-49d6-a72a-0ed3f60d82d6
# ╟─2613110d-1ab8-413a-85ce-a2985ee420da
# ╟─fb53a9ed-df58-410a-8275-e15718514950
# ╠═e5b59cc9-0d14-4d8a-bb25-738540e7ebf9
# ╟─ff27c9fc-0e55-47dc-b189-534f7a48fd3f
# ╟─8ea2593c-2f93-47c1-aa7d-918c848f8bfb
# ╠═69fed6cc-030b-4066-a023-0bbf1637fbbc
# ╠═46cd1ee1-e269-46a7-93d3-72597b53a9a9
# ╠═0e6cf095-5249-4660-8d92-1347b143f795
# ╠═3290dfe4-0ec1-490a-92f9-a2b43d0ae344
# ╠═f3a267f5-f4cd-4056-a16a-538c2ae4d756
# ╠═c04c4d58-0469-45cc-a217-444a2b607245
# ╠═b3f72d4b-9f1f-46fd-8145-212f96c320f8
# ╠═fed9cd8e-2eed-4195-8bab-f6c1b0d3785d
# ╠═1cf7d604-04d9-4b52-a087-aa73c52093e5
# ╟─5e276fd0-887e-4de2-b502-359be36e6fb6
# ╠═4b12f0d9-e4e1-4214-9127-40612f38d7a3
# ╟─9b00810e-8dc8-4602-a185-28e60c027b99
# ╠═8c683b66-1fb2-49ad-9caf-cb891520f5c6
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
# ╠═daa6c126-c5ea-4f9c-b04d-39da54e3fc4c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
