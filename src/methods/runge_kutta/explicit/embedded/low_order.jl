"""
    Fehlberg21(precision::Type{T} = Float64) where T <: AbstractFloat

Fehlberg's second(first)-order method.

https://ntrs.nasa.gov/citations/19690021375
"""
function Fehlberg21(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Fehlberg_2_1
    butcher = SMatrix{4, 5, precision, 20}(
        0, 0, 0, 0,
        1//2, 1//2, 0, 0,
        1, 1//256, 255//256, 0,
        1, 1//512, 255//256, 1//512,
        1, 1//256, 255//256, 0
    ) |> transpose

    ω = SMatrix{2, 3, precision, 6}(
        1, -511//512,
        0, 255//256,
        0, 1//512,
    ) |> transpose

    iteration = Explicit()
    reconstructor = Fehlberg21

    return RungeKutta(name, butcher, iteration, reconstructor; ω)
end

"""
    HeunEuler21(precision::Type{T} = Float64) where T <: AbstractFloat

Heun-Euler second(first)-order method.
"""
function HeunEuler21(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Heun_Euler_2_1
    butcher = SMatrix{3, 4, precision, 12}(
        0, 0, 0,
        1, 1, 0,
        1, 1//2, 1//2,
        1, 1, 0
    ) |> transpose

    ω = SMatrix{2, 2, precision, 4}(
        1, -1//2,
        0, 1//2,
    ) |> transpose

    iteration = Explicit()
    reconstructor = HeunEuler21

    return RungeKutta(name, butcher, iteration, reconstructor; ω)
end

"""
    BogackiShampine32(precision::Type{T} = Float64) where T <: AbstractFloat

Bogacki and Shampine's third(second)-order method.

https://www.sciencedirect.com/science/article/pii/0893965989900797
"""
function BogackiShampine32(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Bogacki_Shampine_3_2
    butcher = SMatrix{5, 6, precision, 30}(
        0, 0, 0, 0, 0,
        1//2, 1//2, 0, 0, 0,
        3//4, 0, 3//4, 0, 0,
        1, 2//9, 1//3, 4//9, 0,
        1, 2//9, 1//3, 4//9, 0,
        1, 7//24, 1//4, 1//3, 1//8
    ) |> transpose

    # note: 3rd order C1 interpolant makes use of FSAL property
    ω = SMatrix{3, 4, precision, 12}(
        1, -4//3, 5//9,
        0, 1, -2//3,
        0, 4//3, -8//9,
        0, -1, 1
    ) |> transpose

    iteration = Explicit()
    reconstructor = BogackiShampine32

    return RungeKutta(name, butcher, iteration, reconstructor; ω)
end
