"""
    AdamsBashforth(; order::Int64, precision::Type{T} = Float64) where T <: AbstractFloat

Adams-Bashforth (AB) explicit multistep method.

Note: `order` ranges from 1-6
"""
function AdamsBashforth(; order::Int64,
                          precision::Type{T} = Float64) where T <: AbstractFloat
    # TODO: tried to find 7th order coefficients online
    @assert 1 <= order <= 6 "order = $order is not valid"

    table = [1 1 0 0 0 0 0
             1 3//2 -1//2 0 0 0 0
             1 23//12 -16//12 5//12 0 0 0
             1 55//24 -59//24 37//24 -9//24 0 0
             1 1901//720 -2774//720 2616//720 -1274//720 251//720 0
             1 4277//1440 -7923//1440 9982//1440 -7298//1440 2877//1440 -475//1440
            ]
    table = table .|> precision
    return Adams(; name = Symbol("Adams_Bashforth_$(order)"), order, table)
end

#=
"""
    AdamsBashforth8(; precision::Type{T} = Float64) where T <: AbstractFloat

Eighth-order Adams-Bashforth multistep method.
"""
function AdamsBashforth8(; precision::Type{T} = Float64) where T <: AbstractFloat
    table = [1 434241//120960 -1152169//120960 2183877//120960 -2664477//120960 2102243//120960 -1041723//120960 295767//120960 -36799//120960]
    table = table .|> precision
    return Adams(; name = :Adams_Bashforth_8, table)
end
=#
