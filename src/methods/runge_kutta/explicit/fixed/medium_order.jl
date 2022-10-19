
# TODO: for Fehlberg4, etc, take Fehlberg45 except 2nd-last (last) row
# TODO: make dictionary with code names
# TODO: make doc strings with links 
# TODO: make documentation page 

function RungeKutta4(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0 0
               1//2 1//2 0 0 0
               1//2 0 1//2 0 0
               1 0 0 1 0
               1 1//6 1//3 1//3 1//6]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Runge_Kutta_4, butcher)
end

function ThreeEightsRule4(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0 0
               1//3 1//3 0 0 0
               2//3 -1//3 1 0 0
               1 1 -1 1 0
               1 1//8 3//8 3//8 1//8]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Three_Eights_Rule_4, butcher)
end

function Ralston4(; precision::Type{<:AbstractFloat} = Float64)
    # TODO: use irrational form
    # note: for multiple terms in element, do a+b+c 0 0 
    # TODO: do sqrt2 = sqrt(precision(2))
    # note: can't use BigFloat precision right now
    if precision == BigFloat
        @warn "Ralston4 can't use BigFloat right now (default to Float64)"
        precision = Float64
    end
    butcher = [0 0 0 0 0
               2//5 2//5 0 0 0
               0.4557372542187894 0.2969776092477536 0.15875964497103556 0 0
               1 0.21810038822592054 -3.050965148692931 3.8328647604670123 0
               1 0.1747602822626904 -0.551480662878733 1.2055355993965235 0.17118478121951902]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Ralston_4, butcher)
end

function Ketcheson4(; precision::Type{<:AbstractFloat} = Float64)
    # SSP
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

    RungeKutta(; name = :Ketcheson_4, butcher)
end

function Butcher5(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0 0 0 0
               1//4 1//4 0 0 0 0 0
               1//4 1//8 1//8 0 0 0 0
               1//2 0 -1//2 1 0 0 0
               3//4 3//16 0 0 9//16 0 0
               1 -3//7 2//7 12//7 -12//7 8//7 0
               1 7//90 0 32//90 12//90 32//90 7//90]
    butcher = butcher .|> precision
     
    RungeKutta(; name = :Butcher_5, butcher)
end

function Butcher6(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0 0 0 0 0
               1//2 1//2 0 0 0 0 0 0
               2//3 2//9 4//9 0 0 0 0 0
               1//3 7//36 2//9 -1//12 0 0 0 0
               5//6 -35//144 -55//36 35//48 15//8 0 0 0
               1//6 -1//360 -11//36 -1//8 1//2 1//10 0 0
               1 -41//260 22//13 43//156 -118//39 32//195 80//39 0
               1 13//200 0 11//40 11//40 4//25 4//25 13//200]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Butcher_6, butcher)
end

