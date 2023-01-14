# TODO: include stability region table/calculator?
"""
Euler's first-order method.
"""
function Euler1(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [0 0
               1 1]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Euler_1, butcher)
end

"""
Heun's second-order method.

Note: strong stability preserving (SSP)
"""
function Heun2(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [0 0 0
               1 1 0
               1 1//2 1//2]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Heun_2, butcher)
end

"""
Second-order midpoint rule.
"""
function Midpoint2(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [0 0 0
               1//2 1//2 0
               1 0 1]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Midpoint_2, butcher)
end

"""
Ralston's second-order method.

https://www.ams.org/journals/mcom/1962-16-080/S0025-5718-1962-0150954-0/S0025-5718-1962-0150954-0.pdf
"""
function Ralston2(; precision::Type{T} = Float64) where T <: AbstractFloat
    # TODO: had load module issues with ChangePrecision
    butcher = [0 0 0
               2//3 2//3 0
               1 1//4 3//4]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Ralston_2, butcher)
end

"""
    Generic2(; alpha::Union{Int, Rational},
               precision::Type{T} = Float64) where T <: AbstractFloat

A generic second-order Runge-Kutta method.

Required parameters: `alpha`
"""
function Generic2(; alpha::Union{Int, Rational},
                    precision::Type{T} = Float64) where T <: AbstractFloat
    @assert alpha != 0 "choose alpha != 0"
    a = alpha

    butcher = [0 0 0
               a a 0
               1 1-1//2a 1//2a]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Generic_2, butcher)
end

"""
Heun's third-order method.
"""
function Heun3(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [0 0 0 0
               1//3 1//3 0 0
               2//3 0 2//3 0
               1 1/4 0 3//4]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Heun_3, butcher)
end

"""
Ralston's third-order method.

https://www.ams.org/journals/mcom/1962-16-080/S0025-5718-1962-0150954-0/S0025-5718-1962-0150954-0.pdf
"""
function Ralston3(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [0 0 0 0
               1//2 1//2 0 0
               3//4 0 3//4 0
               1 2//9 1//3 4//9]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Ralston_3, butcher)
end

"""
Kutta's third-order method.
"""
function RungeKutta3(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [0 0 0 0
               1//2 1//2 0 0
               3//4 0 3//4 0
               1 2//9 1//3 4//9]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Runge_Kutta_3, butcher)
end

"""
Shu and Osher's third-order SSP method.

https://www.sciencedirect.com/science/article/pii/0021999188901775
"""
function ShuOsher3(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [0 0 0 0
               1 1 0 0
               1//2 1//4 1//4 0
               1 1//6 1//6 2//3]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Shu_Osher_3, butcher)
end

"""
Spiteri and Ruuth's third-order SSP method.

https://epubs.siam.org/doi/10.1137/S0036142902419284
"""
function SpiteriRuuth3(; precision::Type{T} = Float64) where T <: AbstractFloat
    # SSP
    butcher = [0 0 0 0 0
               1//2 1//2 0 0 0
               1 1//2 1//2 0 0
               1//2 1//6 1//6 1//6 0
               1 1//6 1//6 1//6 1//2]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Spiteri_Ruuth_3, butcher)
end

"""
    Generic3(; alpha::Union{Int, Rational},
               precision::Type{T} = Float64) where T <: AbstractFloat

A generic third-order Runge-Kutta method.

Required parameters: `alpha`
"""
function Generic3(; alpha::Union{Int, Rational},
                    precision::Type{T} = Float64) where T <: AbstractFloat
    @assert !(alpha in [0, 2//3, 1]) "choose alpha ∉ [0, 2//3, 1]"
    a = alpha

    butcher = [0 0 0 0
               a a 0 0
               1 1+(1-a)//(a*(3a-2)) -(1-a)//(a*(3a-2)) 0
               1 1//2-1//6a 1//(6a*(1-a)) (2-3a)//6(1-a)]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Generic_3, butcher)
end
