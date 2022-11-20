# Linear Algebra

####### Array Initialization ######
# small (column) vector
[1 2]

# small (column) vector
[1 2]

# small array
[1 2; 3 4]

# an array of size (2, 2, 2), not initialized
Vector{Float64}(undef, 2, 2, 2)

# zero array of size (2, 2, 2), the 1st arg is the element type.
zeros(Float64, 2, 2, 2)

# random array
x = randn(Float64, 2, 2, 2)

####### Array Indexing ######
# get the first element.
x[1]

# Yes, it counts from 1!

# get the last element
x[end]
