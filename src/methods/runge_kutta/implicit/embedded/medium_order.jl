
function GaussLegendre42(; precision::Type{T} = Float64) where T <: AbstractFloat
    # note: used BigFloat to reduce float error propagation before .|> precision line
    s = sqrt(BigFloat(3))       # sqrt(3)

    butcher = [1//2-s/6 1//4 1//4-s/6
                1//2+s/6 1//4+s/6 1//4
                1 1//2 1//2
                1 1//2+s/2 1//2-s/2]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Gauss_Legendre_4_2, butcher)
end

function LobattoIIIA42(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [0 0 0 0
               1//2 5//24 1//3 -1//24
               1 1//6 2//3 1//6
               1 1//6 2//3 1//6
               1 -1//2 2 -1//2]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Lobatto_IIIA_4_2, butcher)
end

function LobattoIIIB42(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [0 1//6 -1//6 0
               1//2 1//6 1//3 0
               1 1//6 5//6 0
               1 1//6 2//3 1//6
               1 -1//2 2 -1//2]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Lobatto_IIIB_4_2, butcher)
end

function LobattoIIIC42(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [0 1//6 -1//3 1//6
               1//2 1//6 5//12 -1//12
               1 1//6 2//3 1//6
               1 1//6 2//3 1//6
               1 -1//2 2 -1//2]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Lobatto_IIIC_4_2, butcher)
end

function LobattoIIICS42(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [0 0 0 0
               1//2 1//4 1//4 0
               1 0 1 0
               1 1//6 2//3 1//6
               1 -1//2 2 -1//2]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Lobatto_IIICS_4_2, butcher)
end

function LobattoIIID42(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [0 1//6 0 -1//6
               1//2 1//12 5//12 0
               1 1//2 1//3 1//6
               1 1//6 2//3 1//6
               1 -1//2 2 -1//2]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Lobatto_IIID_4_2, butcher)
end

function RaduaIIA52(; precision::Type{T} = Float64) where T <: AbstractFloat
    # note: used BigFloat to reduce float error propagation before .|> precision line
    s = sqrt(BigFloat(6))       # sqrt(6)

    butcher = [2//5-s/10 11//45-7s/360 37//225-169s/1800 -2//225+s/75
               2//5+s/10 37//225+169s/1800 11//45+7s/360 -2//225-s/75
               1 4//9-s/36 4//9+s/36 1//9
               1 4//9-s/36 4//9+s/36 1//9
               1 1-7s/12 1+7s/12 -1]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Radau_IIA_5_2, butcher)
end

function GaussLegendre64(; precision::Type{T} = Float64) where T <: AbstractFloat
    # note: used BigFloat to reduce float error propagation before .|> precision line
    s = sqrt(BigFloat(15))  # sqrt(15)

    butcher = [1//2-s/10 5//36 2//9-s/15 5//36-s/30
               1//2 5//36+s/24 2//9 5//36-s/24
               1//2+s/10 5//36+s/30 2//9+s/15 5//36
               1 5//18 4//9 5//18
               1 -5//6 8//3 -5//6]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Gauss_Legendre_6_4, butcher)
end
