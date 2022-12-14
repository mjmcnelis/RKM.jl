
function BackwardEuler1(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [1 1
               1 1]
    butcher = butcher .|> precision

    RungeKutta(; name = :Backward_Euler_1, butcher)
end

function ImplicitMidpoint2(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [1//2 1//2
               1 1]
    butcher = butcher .|> precision

    RungeKutta(; name = :Implicit_Midpoint_2, butcher)
end

function QinZhang2(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [1//4 1//4 0
               3//4 1//2 1//4
               1 1//2 1//2]
    butcher = butcher .|> precision

    RungeKutta(; name = :Qin_Zhang_2, butcher)
end

function KraaijevangerSpijker2(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [1//2 1//2 0
               3//2 -1//2 2
               1 -1//2 3//2]
    butcher = butcher .|> precision

    RungeKutta(; name = :Kraaijevanger_Spijker_2, butcher)
end

function PareschiRusso2(; precision::Type{T} = Float64) where T <: AbstractFloat
    # note: used BigFloat to reduce float error propagation before .|> precision line
    s = sqrt(BigFloat(2))   # sqrt(2)
    g = 1 - 1/s             # gamma

    butcher = [g g 0
               1-g 1-2g g
               1 1//2 1//2]
    butcher = butcher .|> precision

    RungeKutta(; name = :Pareschi_Russo_2, butcher)
end

# TODO: maybe use generic formula from
#       https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods
function PareschiRusso3(; precision::Type{T} = Float64) where T <: AbstractFloat
    # note: used BigFloat to reduce float error propagation before .|> precision line
    s = sqrt(BigFloat(2))   # sqrt(2)
    g = 1 - 1/s             # gamma, other options were 1/4 (Qin-Zhang), 1+1/sqrt(2)

    butcher = [g g 0 0
               1-g 1-2g g 0
               1//2 1//2-g 0 g
               1 1//6 1//6 2//3]
    butcher = butcher .|> precision

    RungeKutta(; name = :Pareschi_Russo_3, butcher)
end

function Crouzeix3(; precision::Type{T} = Float64) where T <: AbstractFloat
    # note: used BigFloat to reduce float error propagation before .|> precision line
    s = sqrt(BigFloat(3))   # sqrt(3)

    butcher = [1//2+s/6 1//2+s/6 0
               1//2-s/6 -s/3 1//2+s/6
               1 1//2 1//2]
    butcher = butcher .|> precision

    RungeKutta(; name = :Crouzeix_3, butcher)
end

function RadauIA3(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [0 1//4 -1//4
               2//3 1//4 5//12
               1 1//4 3//4]
    butcher = butcher .|> precision

    RungeKutta(; name = :Radau_IA_3, butcher)
end


function RadauIIA3(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [1//3 5//12 -1//12
               1 3//4 1//4
               1 3//4 1//4]
    butcher = butcher .|> precision

    RungeKutta(; name = :Radau_IIA_3, butcher)
end

# TODO: from https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods
#       but had no specific name
function DIRKL3(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [1//2 1//2 0 0 0
               2//3 1//6 1//2 0 0
               1//2 -1//2 1//2 1//2 0
               1 3//2 -3//2 1//2 1//2
               1 3//2 -3//2 1//2 1//2]
    butcher = butcher .|> precision

    RungeKutta(; name = :DIRK_L_3, butcher)
end
