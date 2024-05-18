
function GaussLegendre42(precision::Type{T} = Float64) where T <: AbstractFloat
    # note: used BigFloat to reduce float error propagation before .|> precision line
    s = sqrt(BigFloat(3))       # sqrt(3)

    name = :Gauss_Legendre_4_2
    butcher = SMatrix{3, 4, precision, 12}(
        1//2-s/6, 1//4, 1//4-s/6,
        1//2+s/6, 1//4+s/6, 1//4,
        1, 1//2, 1//2,
        1, 1//2+s/2, 1//2-s/2
    ) |> transpose
    iteration = FullImplicit()
    reconstructor = GaussLegendre42

    return RungeKutta(name, butcher, iteration, reconstructor)
end

function LobattoIIIA42(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Lobatto_I_I_I_A_4_2
    butcher = SMatrix{4, 5, precision, 20}(
        0, 0, 0, 0,
        1//2, 5//24, 1//3, -1//24,
        1, 1//6, 2//3, 1//6,
        1, 1//6, 2//3, 1//6,
        1, -1//2, 2, -1//2,
    ) |> transpose
    iteration = FullImplicit()
    reconstructor = LobattoIIIA42

    return RungeKutta(name, butcher, iteration, reconstructor)
end

function LobattoIIIB42(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Lobatto_I_I_I_B_4_2
    butcher = SMatrix{4, 5, precision, 20}(
        0, 1//6, -1//6, 0,
        1//2, 1//6, 1//3, 0,
        1, 1//6, 5//6, 0,
        1, 1//6, 2//3, 1//6,
        1, -1//2, 2, -1//2
    ) |> transpose
    iteration = FullImplicit()
    reconstructor = LobattoIIIB42

    return RungeKutta(name, butcher, iteration, reconstructor)
end

function LobattoIIIC42(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Lobatto_I_I_I_C_4_2
    butcher = SMatrix{4, 5, precision, 20}(
        0, 1//6, -1//3, 1//6,
        1//2, 1//6, 5//12, -1//12,
        1, 1//6, 2//3, 1//6,
        1, 1//6, 2//3, 1//6,
        1, -1//2, 2, -1//2
    ) |> transpose
    iteration = FullImplicit()
    reconstructor = LobattoIIIC42

    return RungeKutta(name, butcher, iteration, reconstructor)
end

function LobattoIIICS42(precision::Type{T} = Float64) where T <: AbstractFloat
    # ESDIRK, 3rd stage is explicit, too
    name = :Lobatto_I_I_I_C_S_4_2
    butcher = SMatrix{4, 5, precision, 20}(
        0, 0, 0, 0,
        1//2, 1//4, 1//4, 0,
        1, 0, 1, 0,
        1, 1//6, 2//3, 1//6,
        1, -1//2, 2, -1//2
    ) |> transpose
    iteration = FullImplicit()
    reconstructor = LobattoIIICS42

    return RungeKutta(name, butcher, iteration, reconstructor)
end

function LobattoIIID42(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Lobatto_I_I_I_D_4_2
    butcher = SMatrix{4, 5, precision, 20}(
        0, 1//6, 0, -1//6,
        1//2, 1//12, 5//12, 0,
        1, 1//2, 1//3, 1//6,
        1, 1//6, 2//3, 1//6,
        1, -1//2, 2, -1//2
    ) |> transpose
    iteration = FullImplicit()
    reconstructor = LobattoIIID42

    return RungeKutta(name, butcher, iteration, reconstructor)
end

function RaduaIIA52(precision::Type{T} = Float64) where T <: AbstractFloat
    # note: used BigFloat to reduce float error propagation before .|> precision line
    s = sqrt(BigFloat(6))       # sqrt(6)

    name = :Radau_I_I_A_5_2
    butcher = SMatrix{4, 5, precision, 20}(
        2//5-s/10, 11//45-7s/360, 37//225-169s/1800, -2//225+s/75,
        2//5+s/10, 37//225+169s/1800, 11//45+7s/360, -2//225-s/75,
        1, 4//9-s/36, 4//9+s/36, 1//9,
        1, 4//9-s/36, 4//9+s/36, 1//9,
        1, 1-7s/12, 1+7s/12, -1
    ) |> transpose
    iteration = FullImplicit()
    reconstructor = RaduaIIA52

    return RungeKutta(name, butcher, iteration, reconstructor)
end

function GaussLegendre64(precision::Type{T} = Float64) where T <: AbstractFloat
    # note: used BigFloat to reduce float error propagation before .|> precision line
    s = sqrt(BigFloat(15))  # sqrt(15)

    name = :Gauss_Legendre_6_4
    butcher = SMatrix{4, 5, precision, 20}(
        1//2-s/10, 5//36, 2//9-s/15, 5//36-s/30,
        1//2, 5//36+s/24, 2//9, 5//36-s/24,
        1//2+s/10, 5//36+s/30, 2//9+s/15, 5//36,
        1, 5//18, 4//9, 5//18,
        1, -5//6, 8//3, -5//6
    ) |> transpose
    iteration = FullImplicit()
    reconstructor = GaussLegendre64

    return RungeKutta(name, butcher, iteration, reconstructor)
end
