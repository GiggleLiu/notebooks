#s delay = 5
############  Coding Muscle Training 6: Meta-programming  ############
# Please place your hand on your keyboard, type with me!
# Ready?
# 3. press SPACE to pause.
# 2. press → to move forward.
# 1. press ← to move backward.
# GO!

# functional operations
v = randn(3)

# map(f, v) returns [f(v[1]), f(v[2]), f(v[3])]
map(x->x^2, v)

# reduce(f, v) returns f(f(v[1], v[2]), v[3])
reduce(+, v)

# foldl and foldr are similar to reduce, but with guaranteed left and right associativity.

# mapreduce(r, m, v) returns r(r(m(v[1]), m(v[2])), m(v[3]))
# we can define the l-2 norm function as bellow.
mynorm(v) = sqrt(mapreduce(x->x^2, +, v))

using LinearAlgebra
mynorm(v) ≈ LinearAlgebra.norm(v)

