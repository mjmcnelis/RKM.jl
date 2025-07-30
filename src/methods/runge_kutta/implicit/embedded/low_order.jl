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

"""
    TrapezoidRuleBDF21(precision::Type{T} = Float64) where T <: AbstractFloat

Second-order trapezoid rule BDF method (ABLS-stable) with embedded Euler pair.

R.E. Bank, W.M. Coughran, W. Fichtner, E.H. Grosse, D.J.
Rose, and R.K. Smith. Transient simulation of silicon devices
and circuits. IEEE Transactions on Electron Devices, 32:1992–
2007, 1985.
"""
function TrapezoidRuleBDF21(precision::Type{T} = Float64) where T <: AbstractFloat
    s2 = sqrt(BigFloat(2))   # sqrt(2)
    g = 2 - s2               # gamma

    name = :Trapezoid_Rule_B_D_F_2_1
    butcher = SMatrix{4, 5, precision, 20}(
        0, 0, 0, 0,
        g, g/2, g/2, 0,
        1, 1/(2(2-g)), 1/(2(2-g)), g/2,
        1, 1/(2(2-g)), 1/(2(2-g)), g/2,
        1, 1, 0, 0
    ) |> transpose

    # looks okay (not C1 though)
    ω = SMatrix{2, 3, precision, 6}(
        1 - (1-g)/(4(2-g)),        (3-g)/(4(2-g)) - 1,
        1/(4(2-g)),                1/(4(2-g)),
        -g/(4(2-g)),               1/2 - g/(4(2-g))
    ) |> transpose

    iteration = DiagonalImplicit()
    reconstructor = TrapezoidRuleBDF21

    return RungeKutta(name, butcher, iteration, reconstructor; ω)
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
