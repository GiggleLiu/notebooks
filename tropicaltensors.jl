### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# ╔═╡ 4dd83cba-77ad-11eb-35c1-634ce9675490
using Viznet, Compose

# ╔═╡ b6ac6172-7809-11eb-1df5-c515cba4bd3d
using SimpleTensorNetworks

# ╔═╡ 735d6526-77ad-11eb-2f9c-a9de14c05607
lt = SquareLattice(5, 5);

# ╔═╡ e88cf56a-77af-11eb-0c85-ef3afa7ea1cd
function plotlattice(lt; dx=0.0, dy=0.75)
	Compose.set_default_graphic_size(12cm, 12cm)
	nb = nodestyle(:square)
	eb = bondstyle(:default)
	cb = nodestyle(:circle; r=0.01)
	m, n = size(lt)
	canvas() do
		for i=1:m
			for j=1:n
				j!=n && eb >> (lt[i,j], lt[i,j+1])
				i!=m && eb >> (lt[i,j], lt[i+1,j])
				j!=n && i!=m && eb >> (lt[i+dx,j], lt[i+1-dx,j+1])
				j!=1 && i!=m  && eb >> (lt[i+dx,j], lt[i+1-dx,j-1])
			end
		end
		for i=1:size(lt, 1)
			for j=1:size(lt,2)
				cb >> lt[i,j]
				i!=m && cb >> lt[i+dx,j]
				i!=1 && cb >> lt[i-dx,j]
				i!=m && nb >> lt[i+0.5, j]
				j!=n && nb >> lt[i, j+0.5]
			end
		end
		for i=1:size(lt, 1)
			for j=1:size(lt,2)
				j!=n && i!=m && nb >> lt[i+dy-dx/2, j+dy]
				j!=n && i!=1 && nb >> lt[i-dy+dx/2, j+dy]
			end
		end
	end
end

# ╔═╡ 8108fb92-77b0-11eb-0f7f-b32ef8c0c213
md"The current tensor network mapping..."

# ╔═╡ 31978dfc-77b1-11eb-1fbd-734a8afe0cbd
plotlattice(lt)

# ╔═╡ 6c570198-77b1-11eb-14d5-651ebec2dfbe
md"The correct tensor network should be like this"

# ╔═╡ 5042212c-77b1-11eb-2978-710ffd092a18
plotlattice(lt; dx=0.2)

# ╔═╡ 46586a9c-7809-11eb-1893-a339bdce4000
function δ(::Type{T}, labels...) where T
	data = zeros(T, [2 for l in labels]...)
	data[1] = one(T)
	data[end] = one(T)
	LabeledTensor(data, [labels...])
end

# ╔═╡ a637d580-7809-11eb-30de-9985eed491ca
δ(Float64, 1, 2, 3)

# ╔═╡ 62351f72-780a-11eb-2ba1-376cd7fc9456


# ╔═╡ d238f104-7808-11eb-392c-a5088d7602e7
let n=5
	for i=1:n
		for j=1:n
		end
	end
end

# ╔═╡ Cell order:
# ╠═4dd83cba-77ad-11eb-35c1-634ce9675490
# ╠═735d6526-77ad-11eb-2f9c-a9de14c05607
# ╟─e88cf56a-77af-11eb-0c85-ef3afa7ea1cd
# ╟─8108fb92-77b0-11eb-0f7f-b32ef8c0c213
# ╟─31978dfc-77b1-11eb-1fbd-734a8afe0cbd
# ╟─6c570198-77b1-11eb-14d5-651ebec2dfbe
# ╠═5042212c-77b1-11eb-2978-710ffd092a18
# ╠═b6ac6172-7809-11eb-1df5-c515cba4bd3d
# ╠═46586a9c-7809-11eb-1893-a339bdce4000
# ╠═a637d580-7809-11eb-30de-9985eed491ca
# ╠═62351f72-780a-11eb-2ba1-376cd7fc9456
# ╠═d238f104-7808-11eb-392c-a5088d7602e7
