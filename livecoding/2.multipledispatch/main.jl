######### Function dispatch #########
# dispatch over regular types
f(x::Int) = "integer: $x"
f(3)

# dispatch over type of types
f(x::Type{T}) where T = "type: $T"
f(Mod)

# dispatch over functions
f(x::typeof(f)) = "f"
f(f)