# TODO: for Fehlberg4, etc, take Fehlberg45 except 2nd-last (last) row

"""
    RungeKutta4(precision::Type{T} = Float64) where T <: AbstractFloat

Classic fourth-order Runge-Kutta method.
"""
function RungeKutta4(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Runge_Kutta_4
    butcher = SMatrix{5, 5, precision, 25}(
        0, 0, 0, 0, 0,
        1//2, 1//2, 0, 0, 0,
        1//2, 0, 1//2, 0, 0,
        1, 0, 0, 1, 0,
        1, 1//6, 1//3, 1//3, 1//6
    ) |> transpose
    iteration = Explicit()
    reconstructor = RungeKutta4

    return RungeKutta(name, butcher, iteration, reconstructor)
end

"""
    ThreeEightsRule4(precision::Type{T} = Float64) where T <: AbstractFloat

Fourth-order 3/8 rule.
"""
function ThreeEightsRule4(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Three_Eights_Rule_4
    butcher = SMatrix{5, 5, precision, 25}(
        0, 0, 0, 0, 0,
        1//3, 1//3, 0, 0, 0,
        2//3, -1//3, 1, 0, 0,
        1, 1, -1, 1, 0,
        1, 1//8, 3//8, 3//8, 1//8
    ) |> transpose
    return RungeKutta(name, butcher)
end

"""
    Ralston4(; precision::Type{T} = Float64) where T <: AbstractFloat

Ralston's fourth-order method.

https://www.ams.org/journals/mcom/1962-16-080/S0025-5718-1962-0150954-0/S0025-5718-1962-0150954-0.pdf
"""
function Ralston4(; precision::Type{T} = Float64) where T <: AbstractFloat
    s5 = sqrt(BigFloat(5))       # sqrt(5)

    butcher = [0 0 0 0 0
               2//5 2//5 0 0 0
               7//8-3s5/16 (-2889+1428s5)/1024 (3785-1620s5)/1024 0 0
               1 (-3365+2094s5)/6040 (-975-3046s5)/2552 (467040+203968s5)/240845 0
               1 (263+24s5)/1812 (125-1000s5)/3828 1024(3346+1623s5)/5924787 (30-4s5)/123]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Ralston_4, butcher)
end

"""
    Ketcheson4(; precision::Type{T} = Float64) where T <: AbstractFloat

Ketcheson's fourth-order SSP method.

https://epubs.siam.org/doi/10.1137/07070485X
"""
function Ketcheson4(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [0 0 0 0 0 0 0 0 0 0 0
               1//6 1//6 0 0 0 0 0 0 0 0 0
               1//3 1//6 1//6 0 0 0 0 0 0 0 0
               1//2 1//6 1//6 1//6 0 0 0 0 0 0 0
               2//3 1//6 1//6 1//6 1//6 0 0 0 0 0 0
               1//3 1//15 1//15 1//15 1//15 1//15 0 0 0 0 0
               1//2 1//15 1//15 1//15 1//15 1//15 1//6 0 0 0 0
               2//3 1//15 1//15 1//15 1//15 1//15 1//6 1//6 0 0 0
               5//6 1//15 1//15 1//15 1//15 1//15 1//6 1//6 1//6 0 0
               1 1//15 1//15 1//15 1//15 1//15 1//6 1//6 1//6 1//6 0
               1 1//10 1//10 1//10 1//10 1//10 1//10 1//10 1//10 1//10 1//10]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Ketcheson_4, butcher)
end

"""
    Butcher5(; precision::Type{T} = Float64) where T <: AbstractFloat

Butcher's fifth-order method.
"""
function Butcher5(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [0 0 0 0 0 0 0
               1//4 1//4 0 0 0 0 0
               1//4 1//8 1//8 0 0 0 0
               1//2 0 -1//2 1 0 0 0
               3//4 3//16 0 0 9//16 0 0
               1 -3//7 2//7 12//7 -12//7 8//7 0
               1 7//90 0 32//90 12//90 32//90 7//90]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Butcher_5, butcher)
end

"""
    Butcher6(; precision::Type{T} = Float64) where T <: AbstractFloat

Butcher's sixth-order method.
"""
function Butcher6(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [0 0 0 0 0 0 0 0
               1//2 1//2 0 0 0 0 0 0
               2//3 2//9 4//9 0 0 0 0 0
               1//3 7//36 2//9 -1//12 0 0 0 0
               5//6 -35//144 -55//36 35//48 15//8 0 0 0
               1//6 -1//360 -11//36 -1//8 1//2 1//10 0 0
               1 -41//260 22//13 43//156 -118//39 32//195 80//39 0
               1 13//200 0 11//40 11//40 4//25 4//25 13//200]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Butcher_6, butcher)
end
