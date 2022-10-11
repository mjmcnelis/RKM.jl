
function RungeKutta4(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0 0
               1/2 1/2 0 0 0
               1/2 0 1/2 0 0
               1 0 0 1 0
               1 1/6 1/3 1/3 1/6]
    RungeKutta(:Runge_Kutta_4, butcher .|> precision, Explicit())
end

function ThreeEightsRule4(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0 0
               1/3 1/3 0 0 0
               2/3 -1/3 1 0 0
               1 1 -1 1 0
               1 1/8 3/8 3/8 1/8]
    RungeKutta(:Three_Eights_Rule_4, butcher .|> precision, Explicit())
end

function Ralston4(; precision::Type{<:AbstractFloat} = Float64)
    # TODO: use irrational form
    # note: for multiple terms in element, do a+b+c 0 0 
    butcher = [0 0 0 0 0
               2/5 2/5 0 0 0
               0.4557372542187894 0.2969776092477536 0.15875964497103556 0 0
               1 0.21810038822592054 -3.050965148692931 3.8328647604670123 0
               1 0.1747602822626904 -0.551480662878733 1.2055355993965235 0.17118478121951902]
    RungeKutta(:Ralston_4, butcher .|> precision, Explicit())
end

function Ketcheson4(; precision::Type{<:AbstractFloat} = Float64)
    # SSP
    butcher = [0 0 0 0 0 0 0 0 0 0 0
               1/6 1/6 0 0 0 0 0 0 0 0 0
               1/3 1/6 1/6 0 0 0 0 0 0 0 0
               1/2 1/6 1/6 1/6 0 0 0 0 0 0 0
               2/3 1/6 1/6 1/6 1/6 0 0 0 0 0 0
               1/3 1/15 1/15 1/15 1/15 1/15 0 0 0 0 0
               1/2 1/15 1/15 1/15 1/15 1/15 1/6 0 0 0 0
               2/3 1/15 1/15 1/15 1/15 1/15 1/6 1/6 0 0 0
               5/6 1/15 1/15 1/15 1/15 1/15 1/6 1/6 1/6 0 0
               1 1/15 1/15 1/15 1/15 1/15 1/6 1/6 1/6 1/6 0
               1 1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10 1/10]
    RungeKutta(:Ketcheson_4, butcher .|> precision, Explicit())
end

# TODO: should just grab fehlberg_45 except 2nd last row
fehlberg_4 =[0 0 0 0 0 0;
             1/4 1/4 0 0 0 0;
             3/8 3/32 9/32 0 0 0;
             12/13 1932/2197 -7200/2197 7296/2197 0 0;
             1 439/216 -8 3680/513 -845/4104 0;
             1 25/216 0 1408/2565 2197/4104 -1/5]