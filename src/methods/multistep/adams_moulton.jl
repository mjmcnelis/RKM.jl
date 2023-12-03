"""
    AdamsMoulton1(; precision::Type{T} = Float64) where T <: AbstractFloat

First-order Adams-Moulton multistep method.

Note: second row in `table` are predictor coefficients (i.e. Adams-Bashforth)
"""
function AdamsMoulton1(; precision::Type{T} = Float64) where T <: AbstractFloat
    table = [1 1
             1 1]
    table = table .|> precision
    return LinearMultistep(; name = :Adams_Moulton_1, table)
end

"""
    AdamsMoulton2(; precision::Type{T} = Float64) where T <: AbstractFloat

Second-order Adams-Moulton multistep method.

Note: second row in `table` are predictor coefficients (i.e. Adams-Bashforth)
"""
function AdamsMoulton2(; precision::Type{T} = Float64) where T <: AbstractFloat
    table = [1 1//2 1//2
             1 3//2 -1//2]
    table = table .|> precision
    return LinearMultistep(; name = :Adams_Moulton_2, table)
end

"""
    AdamsMoulton3(; precision::Type{T} = Float64) where T <: AbstractFloat

Third-order Adams-Moulton multistep method.

Note: second row in `table` are predictor coefficients (i.e. Adams-Bashforth)
"""
function AdamsMoulton3(; precision::Type{T} = Float64) where T <: AbstractFloat
    table = [1 5//12 8//12 -1//12
             1 23//12 -16//12 5//12]
    table = table .|> precision
    return LinearMultistep(; name = :Adams_Moulton_3, table)
end

"""
    AdamsMoulton4(; precision::Type{T} = Float64) where T <: AbstractFloat

Fourth-order Adams-Moulton multistep method.

Note: second row in `table` are predictor coefficients (i.e. Adams-Bashforth)
"""
function AdamsMoulton4(; precision::Type{T} = Float64) where T <: AbstractFloat
    table = [1 9//24 19//24 -5//24 1//24
             1 55//24 -59//24 37//24 -9//24]
    table = table .|> precision
    return LinearMultistep(; name = :Adams_Moulton_4, table)
end

"""
    AdamsMoulton5(; precision::Type{T} = Float64) where T <: AbstractFloat

Fifth-order Adams-Moulton multistep method.

Note: second row in `table` are predictor coefficients (i.e. Adams-Bashforth)
"""
function AdamsMoulton5(; precision::Type{T} = Float64) where T <: AbstractFloat
    table = [1 251//720 646//720 -264//720 106//720 -19//720
             1 1901//720 -2774//720 2616//720 -1274//720 251//720]
    table = table .|> precision
    return LinearMultistep(; name = :Adams_Moulton_5, table)
end

"""
    AdamsMoulton6(; precision::Type{T} = Float64) where T <: AbstractFloat

Sixth-order Adams-Moulton multistep method.

Note: second row in `table` are predictor coefficients (i.e. Adams-Bashforth)
"""
function AdamsMoulton6(; precision::Type{T} = Float64) where T <: AbstractFloat
    table = [1 475//1440 1427//1440 -798//1440 482//1440 -173//1440 27//1440
             1 4277//1440 -7923//1440 9982//1440 -7298//1440 2877//1440 -475//1440]
    table = table .|> precision
    return LinearMultistep(; name = :Adams_Moulton_6, table)
end

"""
    AdamsMoulton8(; precision::Type{T} = Float64) where T <: AbstractFloat

Eighth-order Adams-Moulton multistep method.

Note: second row in `table` are predictor coefficients (i.e. Adams-Bashforth)
"""
function AdamsMoulton8(; precision::Type{T} = Float64) where T <: AbstractFloat
    table = [1 36799//120960 139849//120960 -121797//120960 123133//120960 -88547//120960 41499//120960 -11351//120960 1375//120960
             1 434241//120960 -1152169//120960 2183877//120960 -2664477//120960 2102243//120960 -1041723//120960 295767//120960 -36799//120960]
    table = table .|> precision
    return LinearMultistep(; name = :Adams_Moulton_8, table)
end

