#s delay = 5
#########  Coding Muscle Training 1: Basic types and control flow  #########
# Please place your hand on your keyboard, type with me!
# Ready?
# 3. press SPACE to pause.
# 2. press → to move forward.
# 1. press ← to move backward.
# GO!

########### Basic types ###########
# boolean variables
true isa Bool
# size in bytes (1 byte = 8 bits)
sizeof(Bool)
# not
!true
# and
true && false
# or
true || false
# xor: type with \xor<TAB>
true ⊻ false

# integers
3 isa Int
0x3 isa UInt8  # 8 bit, unsigned integer
bitstring(3)
bitstring(0x3)
sizeof(3)
sizeof(0x3)
# operations
7 % 3   # modulus
7 / 3   # / returns a floating point number
7 ÷ 3   # ÷ (by typing: \div<TAB>) returns an integer
7 << 1  # bit shift left
7 >> 1  # bit shift right
7 | 1  # bit-wise or
7 & 1  # bit-wise and
7 ⊻ 1  # bit-wise xor

# floating point numbers and complex numbers
3.0 isa Float64
3f0 isa Float32
3.0 + 3im isa ComplexF64

# operations
3.0 ^ 3   # power

# string tuples and vector
"3.0" isa String
(1, "3.0") isa Tuple
# tuple can be indexed
(1, "3.0")[1]   # the first element

(3.0=>"3.0") isa Pair
# pair has two fields, `first` and `second`
(3.0=>"3.0").first
(3.0=>"3.0").second

# vector
[2, 3.0] isa Vector

# ranges
1:10 isa UnitRange # 1, 2, ..., 10
length(1:10)

1:2:10 isa StepRange

# dict
Dict(3.0=>"3.0") isa Dict

# types can be used for conversion
UInt64(3) isa UInt64
Float32(3) isa Float32
# vectors and tuples can be converted using splatting
[(1, "3.0")...]
([1, 3.0]...,)
((3.0=>"3.0")...,)
(1:2:10...,)
[1:2:10...]
# equivalently
collect(1:2:10)

# the ranges of types
typemin(Int)
typemax(Int)

typemin(UInt64)
typemax(UInt64)

typemin(Float64)
typemax(Float64)

# the zero and one element of types
zero(Float64)
one(Float64)

# What is ∞ * 0?
0.0 * Inf

Inf == Inf

# A special note: NaN is not equal to any value, including itself.
NaN == NaN

# precision of floating point numbers (floating point number distance at 1.0)
eps(Float32)

eps(Float64)

# rounding error is unbiquitus in floating point numbers
1.0 == nextfloat(1.0)

# we prefer using `≈` (typed with \approx<TAB>) in many practical using cases
1.0 ≈ nextfloat(1.0)

########### Arbitary precision ###########
# numbers can overflow
let     # let creates a local scope
    x = 1
    for i=1:100  # iterates over 1, 2, ..., 100
        x *= i
    end
    x   # the (implicitly) returned value
end

# BitInt uses arbitary precision, it avoids overflow
let
    x = BigInt(1)
    for i=1:100
        x *= i
    end
    x
end

########### Control flows ###########
# block statement
# block statement is trivial, it is just a sequence of operations
a = 0; a += 1; a+=2; a

begin
    a = 0
    a += 1
    a+=2
    a
end

# for statement
# the following two statements are the same
for i in 1:3
    println(i)
end

for i = 1:3
    println(i)
end

# declare a variable without initialization
let
    for i = 1:3
        j = i==1 ? 1 : j + 1   # j is not in the outer scope
    end
    j
end

let
    local j   # declared without initialization
    for i = 1:3
        j = i==1 ? 1 : j + 1
    end
    j
end

# using global variables
global_j = 0

let
    for i=1:5
        global_j += 1  # can not access global variables at a local scope
    end
end

let
    for i=1:5
        global global_j += 1  # this one works!
    end
end

global_j


# Special Note: unlike the let statement, the block statement runs in the current scope instead of creating a local one.
#
# if statement
# we can use if statement to find out all prime numbers
using Primes
for i=1:typemax(Int8)
    if isprime(i)
        println(i)
    end
end

# equivalently, you can use a while statement
let i=Int8(0)
    while i < typemax(Int8)
        i += Int8(1)
        if isprime(i)
            println(i)
        end
    end
end

# to find out all prime numbers, you can also use broadcasting 
# first create a range
range = Int8(1):typemax(Int8)
# then create a boolean vector for indexing using broadcasting
boolean_mask = isprime.(range)
# the locations marked with 1 will be collected to a vector
range[boolean_mask]

# even shorter, you can use list comprehension
[x for x in Int8(1):typemax(Int8) if isprime(x)]

# or even shorter with filter
filter(isprime, Int8(1):typemax(Int8))

# check the docstring of filter with the `@doc` macro
@doc filter

# if ... elseif ... else ... end
function compare(x,y)
    if x < y
        relation = "less than"
    elseif x == y
        relation = "equal to"
    else
        relation = "larger than to"
    end
    println("$x is $relation $(y).")  # `$` for interpolation
end

compare(3, 4)

# try ... catch ... finally ... end
try
    x = [1]
    x[0]   # Julia counts from 1!
catch e
    @info "get error $e"  # @info is a fancy way of printing
    throw(e)
finally
    println("I will be executed anyway ;D")
end
