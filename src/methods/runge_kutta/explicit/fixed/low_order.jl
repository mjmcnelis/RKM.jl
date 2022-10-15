# TODO: include stability region table/calculator?
# TODO: BigFloat(1)/BigFloat(3) works but BigFloat(1/3) doesn't 
#       can I rewrite this more elegantly so BigFloat works

function Euler1(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0
               1 1] 
    butcher = butcher .|> precision 

    RungeKutta(; name = :Euler_1, butcher)
end

function Heun2(; precision::Type{<:AbstractFloat} = Float64) 
    # SSP
    butcher = [0 0 0
               1 1 0
               1 1//2 1//2]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Heun_2, butcher)
end

function Midpoint2(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0
               1//2 1//2 0
               1 0 1]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Midpoint_2, butcher)
end

function Ralston2(; precision::Type{<:AbstractFloat} = Float64)
    # TODO: had load module issues with ChangePrecision
    butcher = [0 0 0
               2//3 2//3 0
               1 1//4 3//4]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Ralston_2, butcher)
end

function Generic2(; alpha::Union{Int, Rational}, 
                    precision::Type{<:AbstractFloat} = Float64)
    @assert alpha != 0 "choose alpha != 0"
    a = alpha

    butcher = [0 0 0
               a a 0
               1 1-1//2a 1//2a]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Generic_2, butcher)
end

function Heun3(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0
               1//3 1//3 0 0
               2//3 0 2//3 0
               1 1/4 0 3//4]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Heun_3, butcher)
end

function Ralston3(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0
               1//2 1//2 0 0
               3//4 0 3//4 0
               1 2//9 1//3 4//9]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Ralston_3, butcher)
end

function RungeKutta3(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0
               1//2 1//2 0 0
               3//4 0 3/4 0
               1 2//9 1//3 4//9]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Runge_Kutta_3, butcher)
end

function ShuOsher3(; precision::Type{<:AbstractFloat} = Float64)
    # SSP 
    butcher = [0 0 0 0
               1 1 0 0
               1//2 1//4 1//4 0
               1 1//6 1//6 2//3]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Shu_Osher_3, butcher)
end

function SpiteriRuuth3(; precision::Type{<:AbstractFloat} = Float64)
    # SSP 
    butcher = [0 0 0 0 0
               1//2 1//2 0 0 0
               1 1//2 1//2 0 0
               1//2 1//6 1//6 1//6 0
               1 1//6 1//6 1//6 1//2]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Spiteri_Ruuth_3, butcher)
end

function Generic3(; alpha::Union{Int, Rational}, 
                    precision::Type{<:AbstractFloat} = Float64)
    @assert !(alpha in [0, 2//3, 1]) "choose alpha ∉ [0, 2//3, 1]"
    a = alpha

    butcher = [0 0 0 0 
               a a 0 0 
               1 1+(1-a)//(a*(3a-2)) -(1-a)//(a*(3a-2)) 0
               1 1//2-1//6a 1//(6a*(1-a)) (2-3a)//6(1-a)]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Generic_3, butcher)
end
