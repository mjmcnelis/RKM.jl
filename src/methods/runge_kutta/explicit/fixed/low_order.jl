
# TODO: include stability region table/calculator?

function Euler1(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0
               1 1] 
    RungeKutta(:Euler_1, butcher .|> precision, Explicit())
end

function Heun2(; precision::Type{<:AbstractFloat} = Float64) 
    # SSP
    butcher = [0 0 0
               1 1 0
               1 1/2 1/2]
    RungeKutta(:Heun_2, butcher .|> precision, Explicit())
end

function Midpoint2(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0
               1/2 1/2 0
               1 0 1]
    RungeKutta(:Midpoint_2, butcher .|> precision, Explicit())
end

function Ralston2(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0
               2/3 2/3 0
               1 1/4 3/4]
    RungeKutta(:Ralston_2, butcher .|> precision, Explicit())
end

function Heun3(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0
               1/3 1/3 0 0
               2/3 0 2/3 0
               1 1/4 0 3/4]
    RungeKutta(:Heun_3, butcher .|> precision, Explicit())
end

function Ralston3(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0
               1/2 1/2 0 0
               3/4 0 3/4 0
               1 2/9 1/3 4/9]
    RungeKutta(:Ralston_3, butcher .|> precision, Explicit())
end

function RungeKutta3(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0
               1/2 1/2 0 0
               3/4 0 3/4 0
               1 2/9 1/3 4/9]
    RungeKutta(:Runge_Kutta_3, butcher .|> precision, Explicit())
end

function ShuOsher3(; precision::Type{<:AbstractFloat} = Float64)
    # SSP 
    butcher = [0 0 0 0
               1 1 0 0
               1/2 1/4 1/4 0
               1 1/6 1/6 2/3]
    RungeKutta(:Shu_Osher_3, butcher .|> precision, Explicit())
end

function SpiteriRuuth3(; precision::Type{<:AbstractFloat} = Float64)
    # SSP 
    butcher = [0 0 0 0 0
               1/2 1/2 0 0 0
               1 1/2 1/2 0 0
               1/2 1/6 1/6 1/6 0
               1 1/6 1/6 1/6 1/2]
    RungeKutta(:Spiteri_Ruuth_3, butcher .|> precision, Explicit())
end
