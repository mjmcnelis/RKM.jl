
# TODO: see roots https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods
function Norsett4(; precision::Type{<:AbstractFloat} = Float64)
    s = sqrt(precision(3))      # sqrt(3)
    c = cos(pi/precision(18))   # cos(pi/18)

    x = 1//2 + c/s              # root = 1.068579

    butcher = [x x 0 0
               1//2 1//2-x x 0
               1-x 2x 1-4x x
               1 1/(6(1-2x)^2) 1-1/(3(1-2x)^2) 1/(6(1-2x)^2)]
    butcher = butcher .|> precision

    RungeKutta(; name = :Norsett_4, butcher)
end

function LobattoIIIA4(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0
               1//2 5//24 1//3 -1//24
               1 1//6 2//3 1//6
               1 1//6 2//3 1//6]
    butcher = butcher .|> precision

    RungeKutta(; name = :Lobatto_IIIA_4, butcher)
end

function LobattoIIIB4(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 1//6 -1//6 0
               1//2 1//6 1//3 0
               1 1//6 5//6 0
               1 1//6 2//3 1//6]
    butcher = butcher .|> precision

    RungeKutta(; name = :Lobatto_IIIB_4, butcher)
end

function LobattoIIIC4(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 1//6 -1//3 1//6
               1//2 1//6 5//12 -1//12
               1 1//6 2//3 1//6
               1 1//6 2//3 1//6]
    butcher = butcher .|> precision

    RungeKutta(; name = :Lobatto_IIIC_4, butcher)
end
