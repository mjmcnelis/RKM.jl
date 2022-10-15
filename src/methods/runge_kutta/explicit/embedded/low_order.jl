
# TODO: label Bogacki-Shampine, etc as FSAL

function Fehlberg12(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0
               1//2 1//2 0 0
               1 1//256 255//256 0
               1 1//256 255//256 0
               1 1//512 255//256 1//512]
    EmbeddedRungeKutta(:Fehlberg_12, butcher .|> precision, Explicit())
end

function HeunEuler21(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0
               1 1 0
               1 1//2 1//2
               1 1 0]
    EmbeddedRungeKutta(:Heun_Euler_21, butcher .|> precision, Explicit())
end

function BogackiShampine32(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0 0
               1//2 1//2 0 0 0
               3//4 0 3//4 0 0
               1 2//9 1//3 4//9 0
               1 2//9 1//3 4//9 0
               1 7//24 1//4 1//3 1//8]
    EmbeddedRungeKutta(:Bogacki_Shampine_32, butcher .|> precision, Explicit())
end
