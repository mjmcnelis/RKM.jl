"""
Fehlberg's first(second)-order method.

https://ntrs.nasa.gov/citations/19690021375
"""
function Fehlberg12(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [0 0 0 0
               1//2 1//2 0 0
               1 1//256 255//256 0
               1 1//256 255//256 0
               1 1//512 255//256 1//512]
    butcher = butcher .|> precision

    RungeKutta(; name = :Fehlberg_1_2, butcher)
end

"""
Heun and Euler's second(first)-order method.
"""
function HeunEuler21(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [0 0 0
               1 1 0
               1 1//2 1//2
               1 1 0]
    butcher = butcher .|> precision

    RungeKutta(; name = :Heun_Euler_2_1, butcher)
end

"""
Bogacki and Shampine's third(second)-order method.

https://www.sciencedirect.com/science/article/pii/0893965989900797
"""
function BogackiShampine32(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [0 0 0 0 0
               1//2 1//2 0 0 0
               3//4 0 3//4 0 0
               1 2//9 1//3 4//9 0
               1 2//9 1//3 4//9 0
               1 7//24 1//4 1//3 1//8]
    butcher = butcher .|> precision

    RungeKutta(; name = :Bogacki_Shampine_3_2, butcher)
end
