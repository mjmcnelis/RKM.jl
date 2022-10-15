# TODO: include stability region table/calculator?
# TODO: BigFloat(1)/BigFloat(3) works but BigFloat(1/3) doesn't 
#       can I rewrite this more elegantly so BigFloat works

function Euler1(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0
               1 1] 
    FixedRungeKutta(:Euler_1, butcher .|> precision, Explicit())
end

function Heun2(; precision::Type{<:AbstractFloat} = Float64) 
    # SSP
    butcher = [0 0 0
               1 1 0
               1 1//2 1//2]
    FixedRungeKutta(:Heun_2, butcher .|> precision, Explicit())
end

function Midpoint2(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0
               1//2 1//2 0
               1 0 1]
    FixedRungeKutta(:Midpoint_2, butcher .|> precision, Explicit())
end

function Ralston2(; precision::Type{<:AbstractFloat} = Float64)
    # TODO: had load module issues with ChangePrecision
    butcher = [0 0 0
               2//3 2//3 0
               1 1//4 3//4]
    FixedRungeKutta(:Ralston_2, butcher .|> precision, Explicit())
end

function Generic2(; alpha::Union{Int, Rational}, 
                    precision::Type{<:AbstractFloat} = Float64)
    @assert alpha != 0 "choose alpha != 0"
    a = alpha
    butcher = [0 0 0
               a a 0
               1 1-1//2a 1//2a]
    FixedRungeKutta(:Generic_2, butcher .|> precision, Explicit())
end

function Heun3(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0
               1//3 1//3 0 0
               2//3 0 2//3 0
               1 1/4 0 3//4]
    FixedRungeKutta(:Heun_3, butcher .|> precision, Explicit())
end

function Ralston3(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0
               1//2 1//2 0 0
               3//4 0 3//4 0
               1 2//9 1//3 4//9]
    FixedRungeKutta(:Ralston_3, butcher .|> precision, Explicit())
end

function RungeKutta3(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0
               1//2 1//2 0 0
               3//4 0 3/4 0
               1 2//9 1//3 4//9]
    FixedRungeKutta(:Runge_Kutta_3, butcher .|> precision, Explicit())
end

function ShuOsher3(; precision::Type{<:AbstractFloat} = Float64)
    # SSP 
    butcher = [0 0 0 0
               1 1 0 0
               1//2 1//4 1//4 0
               1 1//6 1//6 2//3]
    FixedRungeKutta(:Shu_Osher_3, butcher .|> precision, Explicit())
end

function SpiteriRuuth3(; precision::Type{<:AbstractFloat} = Float64)
    # SSP 
    butcher = [0 0 0 0 0
               1//2 1//2 0 0 0
               1 1//2 1//2 0 0
               1//2 1//6 1//6 1//6 0
               1 1//6 1//6 1//6 1//2]
    FixedRungeKutta(:Spiteri_Ruuth_3, butcher .|> precision, Explicit())
end

function Generic3(; alpha::Union{Int, Rational}, 
                    precision::Type{<:AbstractFloat} = Float64)
    @assert !(alpha in [0, 2//3, 1]) "choose alpha ∉ [0, 2//3, 1]"
    a = alpha
    butcher = [0 0 0 0 
               a a 0 0 
               1 1+(1-a)//(a*(3a-2)) -(1-a)//(a*(3a-2)) 0
               1 1//2-1//6a 1//(6a*(1-a)) (2-3a)//6(1-a)]
    FixedRungeKutta(:Generic_3, butcher .|> precision, Explicit())
end
