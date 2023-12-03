
function AdamsBashforth1(; precision::Type{T} = Float64) where T <: AbstractFloat
    table = [1 1]
    table = table .|> precision
    return LinearMultistep(; name = :Adams_Bashforth_1, table)
end

function AdamsBashforth2(; precision::Type{T} = Float64) where T <: AbstractFloat
    table = [1 3//2 -1//2]
    table = table .|> precision
    return LinearMultistep(; name = :Adams_Bashforth_2, table)
end

function AdamsBashforth3(; precision::Type{T} = Float64) where T <: AbstractFloat
    table = [1 23//12 -16//12 5//12]
    table = table .|> precision
    return LinearMultistep(; name = :Adams_Bashforth_3, table)
end

function AdamsBashforth4(; precision::Type{T} = Float64) where T <: AbstractFloat
    table = [1 55//24 -59//24 37//24 -9//24]
    table = table .|> precision
    return LinearMultistep(; name = :Adams_Bashforth_4, table)
end

function AdamsBashforth5(; precision::Type{T} = Float64) where T <: AbstractFloat
    table = [1 1901//720 -2774//720 2616//720 -1274//720 251//720]
    table = table .|> precision
    return LinearMultistep(; name = :Adams_Bashforth_5, table)
end

function AdamsBashforth6(; precision::Type{T} = Float64) where T <: AbstractFloat
    table = [1 4277//1440 -7923//1440 9982//1440 -7298//1440 2877//1440 -475//1440]
    table = table .|> precision
    return LinearMultistep(; name = :Adams_Bashforth_6, table)
end

function AdamsBashforth8(; precision::Type{T} = Float64) where T <: AbstractFloat
    table = [1 434241//120960 -1152169//120960 2183877//120960 -2664477//120960 2102243//120960 -1041723//120960 295767//120960 -36799//120960]
    table = table .|> precision
    return LinearMultistep(; name = :Adams_Bashforth_8, table)
end

