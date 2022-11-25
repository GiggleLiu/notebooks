########### Basic types ###########
3 isa Int
3.0 isa Float64
3f0 isa Float32
3.0 + 3im isa ComplexF64
"3.0" isa String
(1, "3.0") isa Tuple
# tuple can be indexed
(1, "3.0")[1]   # the first element

(3.0=>"3.0") isa Pair
# pair has two fields, `first` and `second`
(3.0=>"3.0").first
(3.0=>"3.0").second

Dict(3.0=>"3.0") isa Dict

# types can be used for conversion
Float32(3) isa Float32

########### Functions ###########