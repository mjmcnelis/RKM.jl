# TODO: include stability region table/calculator?
"""
    Euler1(precision::Type{T} = Float64) where T <: AbstractFloat

Euler's first-order method.
"""
function Euler1(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Euler_1
    butcher = SMatrix{2, 2, precision, 4}(
        0, 0,
        1, 1
    ) |> transpose
    iteration = Explicit()
    reconstructor = Euler1

    return RungeKutta(name, butcher, iteration, reconstructor)
end

"""
    Heun2(precision::Type{T} = Float64) where T <: AbstractFloat

Heun's second-order method.

Note: strong stability preserving (SSP)
"""
function Heun2(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Heun_2
    butcher = SMatrix{3, 3, precision, 9}(
        0, 0, 0,
        1, 1, 0,
        1, 1//2, 1//2
    ) |> transpose
    iteration = Explicit()
    reconstructor = Heun2

    return RungeKutta(name, butcher, iteration, reconstructor)
end

"""
    Midpoint2(precision::Type{T} = Float64) where T <: AbstractFloat

Second-order midpoint rule.
"""
function Midpoint2(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Midpoint_2
    butcher = SMatrix{3, 3, precision, 9}(
        0, 0, 0,
        1//2, 1//2, 0,
        1, 0, 1
    ) |> transpose
    iteration = Explicit()
    reconstructor = Midpoint2

    return RungeKutta(name, butcher, iteration, reconstructor)
end

"""
    Ralston2(precision::Type{T} = Float64) where T <: AbstractFloat

Ralston's second-order method.

https://www.ams.org/journals/mcom/1962-16-080/S0025-5718-1962-0150954-0/S0025-5718-1962-0150954-0.pdf
"""
function Ralston2(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Ralston_2
    butcher = SMatrix{3, 3, precision, 9}(
        0, 0, 0,
        2//3, 2//3, 0,
        1, 1//4, 3//4
    ) |> transpose
    iteration = Explicit()
    reconstructor = Ralston2

    return RungeKutta(name, butcher, iteration, reconstructor)
end

"""
    Generic2(precision::Type{T} = Float64;
             alpha::Union{Int, Rational}) where T <: AbstractFloat

A generic second-order Runge-Kutta method.

Required parameters: `alpha`
"""
function Generic2(precision::Type{T} = Float64;
                  alpha::Union{Int, Rational}) where T <: AbstractFloat
    @assert alpha != 0 "choose alpha != 0"
    a = alpha

    # note: reconstruction doesn't work right now (maybe I can make an exception?)
    #       or I could try alpha = rationalize(method.c[2])
    name = :Generic_2
    butcher = SMatrix{3, 3, precision, 9}(
        0, 0, 0,
        a, a, 0,
        1, 1-1//2a, 1//2a
    ) |> transpose
    iteration = Explicit()
    reconstructor = Generic2

    return RungeKutta(name, butcher, iteration, reconstructor)
end

"""
    Heun3(precision::Type{T} = Float64) where T <: AbstractFloat

Heun's third-order method.
"""
function Heun3(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Heun_3
    butcher = SMatrix{4, 4, precision, 16}(
        0, 0, 0, 0,
        1//3, 1//3, 0, 0,
        2//3, 0, 2//3, 0,
        1, 1//4, 0, 3//4
    ) |> transpose
    iteration = Explicit()
    reconstructor = Heun3

    return RungeKutta(name, butcher, iteration, reconstructor)
end

"""
    Ralston3(precision::Type{T} = Float64) where T <: AbstractFloat

Ralston's third-order method.

https://www.ams.org/journals/mcom/1962-16-080/S0025-5718-1962-0150954-0/S0025-5718-1962-0150954-0.pdf
"""
function Ralston3(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Ralston_3
    butcher = SMatrix{4, 4, precision, 16}(
        0, 0, 0, 0,
        1//2, 1//2, 0, 0,
        3//4, 0, 3//4, 0,
        1, 2//9, 1//3, 4//9
    ) |> transpose
    iteration = Explicit()
    reconstructor = Ralston3

    return RungeKutta(name, butcher, iteration, reconstructor)
end

"""
    RungeKutta3(precision::Type{T} = Float64) where T <: AbstractFloat

Kutta's third-order method.
"""
function RungeKutta3(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Runge_Kutta_3
    butcher = SMatrix{4, 4, precision, 16}(
        0, 0, 0, 0,
        1//2, 1//2, 0, 0,
        3//4, 0, 3//4, 0,
        1, 2//9, 1//3, 4//9
    ) |> transpose
    iteration = Explicit()
    reconstructor = RungeKutta3

    return RungeKutta(name, butcher, iteration, reconstructor)
end

"""
    ShuOsher3(precision::Type{T} = Float64) where T <: AbstractFloat

Shu and Osher's third-order SSP method.

https://www.sciencedirect.com/science/article/pii/0021999188901775
"""
function ShuOsher3(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Shu_Osher_3
    butcher = SMatrix{4, 4, precision, 16}(
        0, 0, 0, 0,
        1, 1, 0, 0,
        1//2, 1//4, 1//4, 0,
        1, 1//6, 1//6, 2//3
    ) |> transpose
    iteration = Explicit()
    reconstructor = ShuOsher3

    return RungeKutta(name, butcher, iteration, reconstructor)
end

"""
    SpiteriRuuth3(precision::Type{T} = Float64) where T <: AbstractFloat

Spiteri and Ruuth's third-order SSP method.

https://epubs.siam.org/doi/10.1137/S0036142902419284
"""
function SpiteriRuuth3(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Spiteri_Ruuth_3
    butcher = SMatrix{5, 5, precision, 25}(
        0, 0, 0, 0, 0,
        1//2, 1//2, 0, 0, 0,
        1, 1//2, 1//2, 0, 0,
        1//2, 1//6, 1//6, 1//6, 0,
        1, 1//6, 1//6, 1//6, 1//2
    ) |> transpose
    iteration = Explicit()
    reconstructor = SpiteriRuuth3

    return RungeKutta(name, butcher, iteration, reconstructor)
end

"""
    Generic3(precision::Type{T} = Float64;
             alpha::Union{Int, Rational}) where T <: AbstractFloat

A generic third-order Runge-Kutta method.

Required parameters: `alpha`
"""
function Generic3(precision::Type{T} = Float64;
                  alpha::Union{Int, Rational}) where T <: AbstractFloat
    @assert !(alpha in [0, 2//3, 1]) "choose alpha ∉ [0, 2//3, 1]"
    a = alpha

    butcher = SMatrix{4, 4, precision, 16}(
        0, 0, 0, 0,
        a, a, 0, 0,
        1, 1+(1-a)//(a*(3a-2)), -(1-a)//(a*(3a-2)), 0,
        1, 1//2-1//6a, 1//(6a*(1-a)), (2-3a)//6(1-a)
    ) |> transpose
    iteration = Explicit()
    reconstructor = Generic3

    return RungeKutta(name, butcher, iteration, reconstructor)
end
