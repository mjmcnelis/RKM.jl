"""
    AdamsMoulton(; order::Int64, precision::Type{T} = Float64) where T <: AbstractFloat

Adams-Moulton implicit multistep method.

Note: `order` ranges from 1-6, `table_pred` contains the
       predictor coefficients (i.e. Adams-Bashforth)
"""
function AdamsMoulton(; order::Int64,
                        precision::Type{T} = Float64) where T <: AbstractFloat
    # TODO: tried to find 7th order coefficients online
    @assert 1 <= order <= 6 "order = $order is not valid"

    table = [1 1 0 0 0 0 0
             1 1//2 1//2 0 0 0 0
             1 5//12 8//12 -1//12 0 0 0
             1 9//24 19//24 -5//24 1//24 0 0
             1 251//720 646//720 -264//720 106//720 -19//720 0
             1 475//1440 1427//1440 -798//1440 482//1440 -173//1440 27//1440
            ]
    table = table .|> precision

    table_pred = [1 1 0 0 0 0 0
                  1 3//2 -1//2 0 0 0 0
                  1 23//12 -16//12 5//12 0 0 0
                  1 55//24 -59//24 37//24 -9//24 0 0
                  1 1901//720 -2774//720 2616//720 -1274//720 251//720 0
                  1 4277//1440 -7923//1440 9982//1440 -7298//1440 2877//1440 -475//1440
                ]
    table_pred = table_pred .|> precision

    return Adams(; name = Symbol("Adams_Moulton_$(order)"), order, table, table_pred)
end

#=
"""
    AdamsMoulton8(; precision::Type{T} = Float64) where T <: AbstractFloat

Eighth-order Adams-Moulton multistep method.

Note: second row in `table` are predictor coefficients (i.e. Adams-Bashforth)
"""
function AdamsMoulton8(; precision::Type{T} = Float64) where T <: AbstractFloat
    table = [1 36799//120960 139849//120960 -121797//120960 123133//120960 -88547//120960 41499//120960 -11351//120960 1375//120960
             1 434241//120960 -1152169//120960 2183877//120960 -2664477//120960 2102243//120960 -1041723//120960 295767//120960 -36799//120960]
    table = table .|> precision
    return Adams(; name = :Adams_Moulton_8, table)
end
=#
