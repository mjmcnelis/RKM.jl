"""
    BackwardEuler1(; precision::Type{T} = Float64) where T <: AbstractFloat

First-order backward Euler method (L-stable).
"""
function BackwardEuler1(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [1 1
               1 1]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Backward_Euler_1, butcher)
end

"""
    TrapezoidRuleBDF2(precision::Type{T} = Float64) where T <: AbstractFloat

Second-order trapezoid rule BDF method (ABLS-stable)

R.E. Bank, W.M. Coughran, W. Fichtner, E.H. Grosse, D.J.
Rose, and R.K. Smith. Transient simulation of silicon devices
and circuits. IEEE Transactions on Electron Devices, 32:1992–
2007, 1985.
"""
function TrapezoidRuleBDF2(precision::Type{T} = Float64) where T <: AbstractFloat
    s2 = sqrt(BigFloat(2))   # sqrt(2)
    g = 2 - s2               # gamma
    # note: is it FSAL?
    name = :Trapezoid_Rule_B_D_F_2
    butcher = SMatrix{4, 4, precision, 16}(
        0, 0, 0, 0,
        g, g/2, g/2, 0,
        1, 1/(2(2-g)), 1/(2(2-g)), g/2,
        1, 1/(2(2-g)), 1/(2(2-g)), g/2
    ) |> transpose
    iteration = DiagonalImplicit()
    reconstructor = TrapezoidRuleBDF2

    return RungeKutta(name, butcher, iteration, reconstructor)
end

"""
    ImplicitMidpoint2(; precision::Type{T} = Float64) where T <: AbstractFloat

Second-order implicit mid-point rule (A-stable).
"""
function ImplicitMidpoint2(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [1//2 1//2
               1 1]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Implicit_Midpoint_2, butcher)
end

"""
    QinZhang2(; precision::Type{T} = Float64) where T <: AbstractFloat

Qin and Zhang's econd-order method.
"""
function QinZhang2(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [1//4 1//4 0
               3//4 1//2 1//4
               1 1//2 1//2]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Qin_Zhang_2, butcher)
end

"""
    KraaijevangerSpijker2(; precision::Type{T} = Float64) where T <: AbstractFloat

Kraaijevanger and Spijker's second-order method.
"""
function KraaijevangerSpijker2(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [1//2 1//2 0
               3//2 -1//2 2
               1 -1//2 3//2]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Kraaijevanger_Spijker_2, butcher)
end

"""
    PareschiRusso2(; precision::Type{T} = Float64) where T <: AbstractFloat

Pareschi and Russo's second-order method.
"""
function PareschiRusso2(; precision::Type{T} = Float64) where T <: AbstractFloat
    s2 = sqrt(BigFloat(2))   # sqrt(2)
    g = 1 - 1/s2             # gamma

    butcher = [g g 0
               1-g 1-2g g
               1 1//2 1//2]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Pareschi_Russo_2, butcher)
end

"""
    PareschiRusso3(; precision::Type{T} = Float64) where T <: AbstractFloat

Pareschi and Russo's third-order method.
"""
function PareschiRusso3(; precision::Type{T} = Float64) where T <: AbstractFloat
    # TODO: maybe use generic formula from
    #       https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods
    s2 = sqrt(BigFloat(2))   # sqrt(2)
    g = 1 - 1/s2             # gamma, other options were 1/4 (Qin-Zhang), 1+1/sqrt(2)

    butcher = [g g 0 0
               1-g 1-2g g 0
               1//2 1//2-g 0 g
               1 1//6 1//6 2//3]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Pareschi_Russo_3, butcher)
end

"""
    Crouzeix3(; precision::Type{T} = Float64) where T <: AbstractFloat

Crouzeix's third-order method.
"""
function Crouzeix3(; precision::Type{T} = Float64) where T <: AbstractFloat
    s3 = sqrt(BigFloat(3))   # sqrt(3)

    butcher = [1//2+s3/6 1//2+s3/6 0
               1//2-s3/6 -s3/3 1//2+s3/6
               1 1//2 1//2]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Crouzeix_3, butcher)
end

"""
    RadauIA3(; precision::Type{T} = Float64) where T <: AbstractFloat

Radau IA3 third-order method.
"""
function RadauIA3(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [0 1//4 -1//4
               2//3 1//4 5//12
               1 1//4 3//4]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Radau_IA_3, butcher)
end

"""
    RadauIIA3(; precision::Type{T} = Float64) where T <: AbstractFloat

Radau IIA3 third-order method.
"""
function RadauIIA3(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [1//3 5//12 -1//12
               1 3//4 1//4
               1 3//4 1//4]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Radau_IIA_3, butcher)
end


"""
    DIRKL3(; precision::Type{T} = Float64) where T <: AbstractFloat

Third-order L-stable diagonal implicit method.

https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods
"""
function DIRKL3(; precision::Type{T} = Float64) where T <: AbstractFloat

    # note: method is FSAL but not first-stage explicit, but it's possible
    #       I could use FSAL to initialize guess (not  very high priority)
    # TODO: find more specific name
    butcher = [1//2 1//2 0 0 0
               2//3 1//6 1//2 0 0
               1//2 -1//2 1//2 1//2 0
               1 3//2 -3//2 1//2 1//2
               1 3//2 -3//2 1//2 1//2]
    butcher = butcher .|> precision

    return RungeKutta(; name = :DIRK_L_3, butcher)
end
