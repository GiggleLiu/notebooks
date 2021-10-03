# https://susan-stepney.blogspot.com/2014/02/mathjax.html
using Yao
using LuxurySparse: PermMatrix, IMatrix
using Expokit
using Plots

function quantumwalk(N::Int, t::Real)
    P = 2N + 1
    T = ComplexF64
    C = kron(mat(T, H), IMatrix{P}())
    shiftleft = PermMatrix(circshift(1:P, -1), ones(T, P))
    shiftright = PermMatrix(circshift(1:P, 1), ones(T, P))
    S = kron(mat(T, ConstGate.P0), shiftleft) + kron(mat(T, ConstGate.P1), shiftright)
    # symmetric initial state, use `im` to avoid left-right movers interfere.
    reg = kron([1.0, 1.0im]/sqrt(2), (x=zeros(T, P); x[N+1]=1; x))

    for i=1:t
        reg = S * C * reg
    end

    plot(sum(abs2.(reshape(reg, :, 2)), dims=2) |> vec)
end

quantumwalk(100, 150)

function quantumwalk2(N::Int, t::Real)
    P = 2N + 1
    T = ComplexF64
    shiftleft = PermMatrix(circshift(1:P, -1), ones(T, P))
    shiftright = PermMatrix(circshift(1:P, 1), ones(T, P))
    H = shiftleft + shiftright
    # symmetric initial state, use `im` to avoid left-right movers interfere.
    reg = zeros(T, P); reg[N+1]=1
    reg = Expokit.expmv(-im*t, H, reg)
    plot(abs2.(reg) |> vec)
end

quantumwalk2(100, 40)