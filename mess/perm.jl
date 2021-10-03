using Base.Cartesian
using Base: size_to_strides, checkdims_perm
using Random

function myperm!(P, B, perm)
    checkdims_perm(P, B, perm)
    native_strides = size_to_strides(1, size(B)...)
    strides_1 = 0
    begin
        strides_2 = native_strides[perm[1]]
        strides_3 = native_strides[perm[2]]
        strides_4 = native_strides[perm[3]]
        strides_5 = native_strides[perm[4]]
        strides_6 = native_strides[perm[5]]
        strides_7 = native_strides[perm[6]]
        strides_8 = native_strides[perm[7]]
        strides_9 = native_strides[perm[8]]
        strides_10 = native_strides[perm[9]]
        strides_11 = native_strides[perm[10]]
        strides_12 = native_strides[perm[11]]
        strides_13 = native_strides[perm[12]]
        strides_14 = native_strides[perm[13]]
        strides_15 = native_strides[perm[14]]
        strides_16 = native_strides[perm[15]]
        strides_17 = native_strides[perm[16]]
        strides_18 = native_strides[perm[17]]
        strides_19 = native_strides[perm[18]]
        strides_20 = native_strides[perm[19]]
        strides_21 = native_strides[perm[20]]
        strides_22 = native_strides[perm[21]]
        strides_23 = native_strides[perm[22]]
        strides_24 = native_strides[perm[23]]
        strides_25 = native_strides[perm[24]]
        strides_26 = native_strides[perm[25]]
    end
    offset = 1 - sum((strides_2, strides_3, strides_4, strides_5, strides_6, strides_7, strides_8, strides_9, strides_10, strides_11, strides_12, strides_13, strides_14, strides_15, strides_16, strides_17, strides_18, strides_19, strides_20, strides_21, strides_22, strides_23, strides_24, strides_25, strides_26))
    ind = 1
    begin
        counts_26 = strides_26
    end
    begin
        for i_25 = Base.axes(P, 25)
            counts_25 = strides_25
            begin
                for i_24 = Base.axes(P, 24)
                    counts_24 = strides_24
                    begin
                        for i_23 = Base.axes(P, 23)
                            counts_23 = strides_23
                            begin
                                for i_22 = Base.axes(P, 22)
                                    counts_22 = strides_22
                                    begin
                                        for i_21 = Base.axes(P, 21)
                                            counts_21 = strides_21
                                            begin
                                                for i_20 = Base.axes(P, 20)
                                                    counts_20 = strides_20
                                                    begin
                                                        for i_19 = Base.axes(P, 19)
                                                            counts_19 = strides_19
                                                            begin
                                                                for i_18 = Base.axes(P, 18)
                                                                    counts_18 = strides_18
                                                                    begin
                                                                        for i_17 = Base.axes(P, 17)
                                                                            counts_17 = strides_17
                                                                            begin
                                                                                for i_16 = Base.axes(P, 16)
                                                                                    counts_16 = strides_16
                                                                                    begin
                                                                                        for i_15 = Base.axes(P, 15)
                                                                                            counts_15 = strides_15
                                                                                            begin
                                                                                                for i_14 = Base.axes(P, 14)
                                                                                                    counts_14 = strides_14
                                                                                                    begin
                                                                                                        for i_13 = Base.axes(P, 13)
                                                                                                            counts_13 = strides_13
                                                                                                            begin
                                                                                                                for i_12 = Base.axes(P, 12)
                                                                                                                    counts_12 = strides_12
                                                                                                                    begin
                                                                                                                        for i_11 = Base.axes(P, 11)
                                                                                                                            counts_11 = strides_11
                                                                                                                            begin
                                                                                                                                for i_10 = Base.axes(P, 10)
                                                                                                                                    counts_10 = strides_10
                                                                                                                                    begin
                                                                                                                                        for i_9 = Base.axes(P, 9)
                                                                                                                                            counts_9 = strides_9
                                                                                                                                            begin
                                                                                                                                                for i_8 = Base.axes(P, 8)
                                                                                                                                                    counts_8 = strides_8
                                                                                                                                                    begin
                                                                                                                                                        for i_7 = Base.axes(P, 7)
                                                                                                                                                            counts_7 = strides_7
                                                                                                                                                            begin
                                                                                                                                                                for i_6 = Base.axes(P, 6)
                                                                                                                                                                    counts_6 = strides_6
                                                                                                                                                                    begin
                                                                                                                                                                        for i_5 = Base.axes(P, 5)
                                                                                                                                                                            counts_5 = strides_5
                                                                                                                                                                            begin
                                                                                                                                                                                for i_4 = Base.axes(P, 4)
                                                                                                                                                                                    counts_4 = strides_4
                                                                                                                                                                                    begin
                                                                                                                                                                                        for i_3 = Base.axes(P, 3)
                                                                                                                                                                                            counts_3 = strides_3
                                                                                                                                                                                            begin
                                                                                                                                                                                                for i_2 = Base.axes(P, 2)
                                                                                                                                                                                                    counts_2 = strides_2
                                                                                                                                                                                                    begin
                                                                                                                                                                                                        for i_1 = Base.axes(P, 1)
                                                                                                                                                                                                            counts_1 = strides_1
                                                                                                                                                                                                            begin
                                                                                                                                                                                                                sumc = sum([counts_2, counts_3, counts_4, counts_5, counts_6, counts_7, counts_8, counts_9, counts_10, counts_11, counts_12, counts_13, counts_14, counts_15, counts_16, counts_17, counts_18, counts_19, counts_20, counts_21, counts_22, counts_23, counts_24, counts_25, counts_26])
                                                                                                                                                                                                                @inbounds P[ind] = B[sumc + offset]
                                                                                                                                                                                                                ind += 1
                                                                                                                                                                                                            end
                                                                                                                                                                                                            counts_2 += strides_2
                                                                                                                                                                                                        end
                                                                                                                                                                                                    end
                                                                                                                                                                                                    counts_3 += strides_3
                                                                                                                                                                                                end
                                                                                                                                                                                            end
                                                                                                                                                                                            counts_4 += strides_4
                                                                                                                                                                                        end
                                                                                                                                                                                    end
                                                                                                                                                                                    counts_5 += strides_5
                                                                                                                                                                                end
                                                                                                                                                                            end
                                                                                                                                                                            counts_6 += strides_6
                                                                                                                                                                        end
                                                                                                                                                                    end
                                                                                                                                                                    counts_7 += strides_7
                                                                                                                                                                end
                                                                                                                                                            end
                                                                                                                                                            counts_8 += strides_8
                                                                                                                                                        end
                                                                                                                                                    end
                                                                                                                                                    counts_9 += strides_9
                                                                                                                                                end
                                                                                                                                            end
                                                                                                                                            counts_10 += strides_10
                                                                                                                                        end
                                                                                                                                    end
                                                                                                                                    counts_11 += strides_11
                                                                                                                                end
                                                                                                                            end
                                                                                                                            counts_12 += strides_12
                                                                                                                        end
                                                                                                                    end
                                                                                                                    counts_13 += strides_13
                                                                                                                end
                                                                                                            end
                                                                                                            counts_14 += strides_14
                                                                                                        end
                                                                                                    end
                                                                                                    counts_15 += strides_15
                                                                                                end
                                                                                            end
                                                                                            counts_16 += strides_16
                                                                                        end
                                                                                    end
                                                                                    counts_17 += strides_17
                                                                                end
                                                                            end
                                                                            counts_18 += strides_18
                                                                        end
                                                                    end
                                                                    counts_19 += strides_19
                                                                end
                                                            end
                                                            counts_20 += strides_20
                                                        end
                                                    end
                                                    counts_21 += strides_21
                                                end
                                            end
                                            counts_22 += strides_22
                                        end
                                    end
                                    counts_23 += strides_23
                                end
                            end
                            counts_24 += strides_24
                        end
                    end
                    counts_25 += strides_25
                end
            end
            counts_26 += strides_26
        end
    end
    return P
end

#(n=25; p=randn(fill(2, n)...); t = randn(fill(2, n)...); @time myperm!(p, t, randperm(n)));


for (V, PT, BT) in Any[((:N,), BitArray, BitArray), ((:T,:N), Array, StridedArray)]
    @eval @generated function newperm!(P::$PT{$(V...)}, B::$BT{$(V...)}, perm) where $(V...)
        quote
            checkdims_perm(P, B, perm)

            #calculates all the strides
            native_strides = size_to_strides(1, size(B)...)
            strides_1 = 0
            @nexprs $N d->(strides_{d+1} = native_strides[perm[d]])

            #Creates offset, because indexing starts at 1
            offset = 1 - sum(@ntuple $N d->strides_{d+1})

            sumc = 0
            ind = 1
            @nexprs 1 d->(counts_{$N+1} = strides_{$N+1}) # a trick to set counts_($N+1)
            @nloops($N, i, P,
                    d->(df_d=i_d*strides_{d+1} ;sumc += df_d), # PRE
                    d->(sumc -= df_d), # POST
                    begin # BODY
                        @inbounds P[ind] = B[sumc+offset]
                        ind += 1
                    end)

            return P
        end
    end
end


using Test

@testset "newperm" begin
    n=25
    t=randn(rand(1:2, n)...)
    perm = randperm(n)
    p = zeros(eltype(t), size.(Ref(t), (perm...,)));
    @time newperm!(p, t, perm);
    # 0.395072 seconds (520.17 k allocations: 20.894 MiB, 99.99% compilation time)
    @time permutedims!(p, t, perm);
    # 41.520901 seconds (502.11 k allocations: 20.155 MiB, 100.00% compilation time)
    @test newperm!(p, t, perm) ≈ permutedims!(p, t, perm)
end

using BenchmarkTools
# high dim
n=25
t=randn(rand(2:2, n)...)
perm = randperm(n)
p = zeros(eltype(t), size.(Ref(t), (perm...,)))
@btime newperm!($p, $t, $perm)
# 179.814 ms (2 allocations: 96 bytes)
@btime permutedims!($p, $t, $perm)
# 201.038 ms (2 allocations: 96 bytes)

# low dim (make sure no performance regression)
t=randn(100, 200, 300)
perm = [3,2,1]
p = zeros(eltype(t), size.(Ref(t), (perm...,)))
@btime newperm!($p, $t, $perm)
# 19.203 ms (2 allocations: 96 bytes)
@btime permutedims!($p, $t, $perm);
# 19.339 ms (2 allocations: 96 bytes)