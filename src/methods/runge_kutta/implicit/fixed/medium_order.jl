
# TODO: see roots https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods
function Norsett4(; precision::Type{T} = Float64) where T <: AbstractFloat
    # note: used BigFloat to reduce float error propagation before .|> precision line
    s = sqrt(BigFloat(3))       # sqrt(3)
    c = cos(pi/BigFloat(18))    # cos(pi/18)

    x = 1//2 + c/s              # root = 1.068579 was recommended

    butcher = [x x 0 0
               1//2 1//2-x x 0
               1-x 2x 1-4x x
               1 1/(6(1-2x)^2) 1-1/(3(1-2x)^2) 1/(6(1-2x)^2)]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Norsett_4, butcher)
end

function RaduaIA5(; precision::Type{T} = Float64) where T <: AbstractFloat
    # note: used BigFloat to reduce float error propagation before .|> precision line
    s = sqrt(BigFloat(6))       # sqrt(6)

    butcher = [0 1//9 -(1+s)/18 (-1+s)/18
               3//5-s/10 1//9 11//45+7s/360 11//45-43s/360
               3//5+s/10 1//9 11//45+43s/360 11//45-7s/360
               1 1//9 4//9+s/36 4//9-s/36]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Radau_IA_5, butcher)
end
