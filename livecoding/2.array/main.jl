#s delay = 5
############  Coding Muscle Training 2: Array Operations  ############
# Please place your hand on your keyboard, type with me!
# Ready?
# 3
# 2
# 1
# GO!

####### Part 0: Array Types ######
# All array types are subtypes of the abstract type `AbstractArray`
# Let us print its type-tree
using AbstractTrees

AbstractTrees.children(x::Type) = subtypes(x)

print_tree(AbstractArray)

# If you dump the concrete type `Array`
dump(Array)

# You will see this type has two type parameters
# The first one is the element type
eltype(Array{Float64, 4})

# The second one is the tensor rank
ndims(Array{Float64, 4})

####### Part 1: Vectors ######
# initialized to zeros
v_zero = zeros(Float64, 3)

# Vector is an alias of `Array` type
typeof(v_zero)

# initialized to undefined values
v_undef = Vector{Float64}(undef, 3)

# initialized to uniform random values
v_rand = rand(Float64, 3)  # the first argument is optional

# initialized to normal random values
v_randn = randn(3)

# initialized to 1..n
v_range = collect(1:3)

# which is equivalent to
for i=1:3
    v_undef[i] = i  # setting value
end
v_undef

# manual specification
v = [1, 2, 3]

# length of a vector
length(v)

# get the first element, Julia arrays counts from 1!

v[1]

# get the last element
v[end]

# add an element at the end
push!(v, 4)

# insert an element at location 1, with value 0
insert!(v, 1, 0)

# collect elements at location 2, 3 and 4 into a new vector (allocates!)
v_allocate = v[2:4]

# v and v_allocate do not share memory
v_allocate[2] = 10

v

# get a subarray without allocation
v_view = view(v, 2:4)

# v and v_view share memory
v_view[2] = 10

v

# column vector with content type Any
v_any = Any[1, 2, 3]

# operating a generic typed vector can be much slower.
using BenchmarkTools
@btime sum($v_any)

@btime sum($v)

# you cannot add a string element into a vector of non-string element type
push!(v, "muscle")

# but it is totally fine for a vector with generic type
push!(v_any, "muscle")

####### Part 2: Matrices and Tensors ######
# an tensor of size (2, 2, 2) wtih uninitialized content
Array{Float64}(undef, 2, 2, 2)

# an tensor of size (2, 2, 2) innitialized to all zero, the 1st arg is the element type.
zeros(Float64, 2, 2, 2)

# matrix isa a alias of Array{T, 2}
m = [1 2; 3 4]

# get the size of a matrix
size(m)

# get the first row
m[1:1, :]

# tensors can also be indexed like vectors
m[1]

# one should be extremely careful that Julia arrays are indexed in a column major order (unlike C and python).
m[2]

# matrix multiplication
m * m

# broadcasted multiplication
m .* m

# broadcasted operations are faster than doing operations separately.
M = randn(100, 100);

@btime M .^ 2 .* 2;

@btime (M .^ 2) * 2;

# this is because broadcasting fuses loops and avoids some allocations.

# row vector
v_row = [1 2]

# it is equivalent to transposing a column vector
v_col = [1, 2]
v_row == v_col'

# vector matrix multiplication
v_row * m

# for complex number, ' means hermitian conjugate.
m = randn(ComplexF64, 2, 2)

m'