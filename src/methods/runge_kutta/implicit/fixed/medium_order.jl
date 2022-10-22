
# TODO: see roots https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods
function Norsett4(; precision::Type{<:AbstractFloat} = Float64)
    # note: used BigFloat to reduce float error propagation before .|> precision line
    s = sqrt(BigFloat(3))       # sqrt(3)
    c = cos(pi/BigFloat(18))    # cos(pi/18)

    x = 1//2 + c/s              # root = 1.068579 was recommended

    butcher = [x x 0 0
               1//2 1//2-x x 0
               1-x 2x 1-4x x
               1 1/(6(1-2x)^2) 1-1/(3(1-2x)^2) 1/(6(1-2x)^2)]
    butcher = butcher .|> precision

    RungeKutta(; name = :Norsett_4, butcher)
end

function GaussLegendre4(; precision::Type{<:AbstractFloat} = Float64)
    # note: used BigFloat to reduce float error propagation before .|> precision line
    s = sqrt(BigFloat(3))       # sqrt(3)

    butcher = [1//2-s/6 1//4 1//4-s/6
               1//2+s/6 1//4+s/6 1//4
               1 1//2 1//2]
    butcher = butcher .|> precision

    RungeKutta(; name = :Gauss_Legendre_4, butcher)
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

function LobattoIIICS4(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0
               1//2 1//4 1//4 0
               1 0 1 0
               1 1//6 2//3 1//6]
    butcher = butcher .|> precision

    RungeKutta(; name = :Lobatto_IIICS_4, butcher)
end

function LobattoIIID4(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 1//6 0 -1//6
               1//2 1//12 5//12 0
               1 1//2 1//3 1//6
               1 1//6 2//3 1//6]
    butcher = butcher .|> precision

    RungeKutta(; name = :Lobatto_IIID_4, butcher)
end

function RaduaIA5(; precision::Type{<:AbstractFloat} = Float64)
    # note: used BigFloat to reduce float error propagation before .|> precision line
    s = sqrt(BigFloat(6))       # sqrt(6)

    butcher = [0 1//9 -(1+s)/18 (-1+s)/18
               3//5-s/10 1//9 11//45+7s/360 11//45-43s/360
               3//5+s/10 1//9 11//45+43s/360 11//45-7s/360
               1 1//9 4//9+s/36 4//9-s/36]
    butcher = butcher .|> precision

    RungeKutta(; name = :Radau_IA_5, butcher)
end

function RaduaIIA5(; precision::Type{<:AbstractFloat} = Float64)
    # note: used BigFloat to reduce float error propagation before .|> precision line
    s = sqrt(BigFloat(6))       # sqrt(6)

    butcher = [2//5-s/10 11//45-7s/360 37//225-169s/1800 -2//225+s/75
               2//5+s/10 37//225+169s/1800 11//45+7s/360 -2//225-s/75
               1 4//9-s/36 4//9+s/36 1//9
               1 4//9-s/36 4//9+s/36 1//9]
    butcher = butcher .|> precision

    RungeKutta(; name = :Radau_IIA_5, butcher)
end

function GaussLegendre6(; precision::Type{<:AbstractFloat} = Float64)
    # note: used BigFloat to reduce float error propagation before .|> precision line
    s = sqrt(BigFloat(15))  # sqrt(15)

    butcher = [1//2-s/10 5//36 2//9-s/15 5//36-s/30
               1//2 5//36+s/24 2//9 5//36-s/24
               1//2+s/10 5//36+s/30 2//9+s/15 5//36
               1 5//18 4//9 5//18]
    butcher = butcher .|> precision

    RungeKutta(; name = :Gauss_Legendre_6, butcher)
end
