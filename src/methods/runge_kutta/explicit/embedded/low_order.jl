
function Fehlberg12(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0
               1//2 1//2 0 0
               1 1//256 255//256 0
               1 1//256 255//256 0
               1 1//512 255//256 1//512]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Fehlberg_12, butcher)
end

function HeunEuler21(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0
               1 1 0
               1 1//2 1//2
               1 1 0]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Heun_Euler_21, butcher)
end

function BogackiShampine32(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0 0
               1//2 1//2 0 0 0
               3//4 0 3//4 0 0
               1 2//9 1//3 4//9 0
               1 2//9 1//3 4//9 0
               1 7//24 1//4 1//3 1//8]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Bogacki_Shampine_32, butcher)
end
