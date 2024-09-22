"""
    ImplicitTrapezoid21(precision::Type{T} = Float64) where T <: AbstractFloat

Second-order implicit trapezoid rule (A-stable) with embedded Euler pair.
"""
function ImplicitTrapezoid21(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Implicit_Trapezoid_2_1
    butcher = SMatrix{3, 4, precision, 12}(
        0, 0, 0,
        1, 1//2, 1//2,
        1, 1//2, 1//2,
        1, 1, 0
    ) |> transpose
    iteration = DiagonalImplicit()
    reconstructor = ImplicitTrapezoid21

    return RungeKutta(name, butcher, iteration, reconstructor)
end

function LobattoIIIB21(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Lobatto_I_I_I_B_2_1
    butcher = SMatrix{3, 4, precision, 12}(
        1//2, 1//2, 0,
        1//2, 1//2, 0,
        1, 1//2, 1//2,
        1, 1, 0
    ) |> transpose
    iteration = DiagonalImplicit()
    reconstructor = LobattoIIIB21

    return RungeKutta(name, butcher, iteration, reconstructor)
end

function LobattoIIIC21(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Lobatto_I_I_I_C_2_1
    butcher = SMatrix{3, 4, precision, 12}(
        0, 1//2, -1//2,
        1, 1//2, 1//2,
        1, 1//2, 1//2,
        1, 1, 0
    ) |> transpose
    iteration = FullImplicit()
    reconstructor = LobattoIIIC21

    return RungeKutta(name, butcher, iteration, reconstructor)
end
