"""
    Euler1(precision::Type{T} = Float64) where T <: AbstractFloat

Euler's explicit first-order method.
"""
function Euler1(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Euler_1
    butcher = SMatrix{2, 2, precision, 4}(
        0, 0,
        1, 1
    ) |> transpose

    ω = SMatrix{1, 1, precision, 1}(
        1
    ) |> transpose

    iteration = Explicit()
    reconstructor = Euler1

    return RungeKutta(name, butcher, iteration, reconstructor; ω)
end

"""
    Heun2(precision::Type{T} = Float64) where T <: AbstractFloat

Heun's second-order method.

Note: strong stability preserving (SSP)
"""
function Heun2(precision::Type{T} = Float64) where T <: AbstractFloat
    # TODO: don't use name for reconstruction anymore
    #       so should set order = [2.0, 1.0] directly
    name = :Heun_2_1
    butcher = SMatrix{3, 4, precision, 12}(
        0, 0, 0,
        1, 1, 0,
        1, 1//2, 1//2,
        1, 1, 0
    ) |> transpose

    ω = SMatrix{2, 2, precision, 4}(
        1, -1//2,
        0, 1//2
    ) |> transpose

    iteration = Explicit()
    reconstructor = Heun2

    return RungeKutta(name, butcher, iteration, reconstructor; ω)
end

"""
    Midpoint2(precision::Type{T} = Float64) where T <: AbstractFloat

Second-order midpoint rule.
"""
function Midpoint2(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Midpoint_2_1
    butcher = SMatrix{3, 4, precision, 12}(
        0, 0, 0,
        1//2, 1//2, 0,
        1, 0, 1,
        1, 1, 0
    ) |> transpose

    ω = SMatrix{2, 2, precision, 4}(
        1, -1,
        0, 1
    ) |> transpose

    iteration = Explicit()
    reconstructor = Midpoint2

    return RungeKutta(name, butcher, iteration, reconstructor; ω)
end

"""
    Ralston2(precision::Type{T} = Float64) where T <: AbstractFloat

Ralston's second-order method.

https://www.ams.org/journals/mcom/1962-16-080/S0025-5718-1962-0150954-0/S0025-5718-1962-0150954-0.pdf
"""
function Ralston2(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Ralston_2_1
    butcher = SMatrix{3, 4, precision, 12}(
        0, 0, 0,
        2//3, 2//3, 0,
        1, 1//4, 3//4,
        1, 1, 0
    ) |> transpose

    ω = SMatrix{2, 2, precision, 4}(
        1, -3//4,
        0, 3//4
    ) |> transpose

    iteration = Explicit()
    reconstructor = Ralston2

    return RungeKutta(name, butcher, iteration, reconstructor; ω)
end

"""
    Fehlberg2(precision::Type{T} = Float64) where T <: AbstractFloat

Fehlberg's second-order method.

https://ntrs.nasa.gov/citations/19690021375
"""
function Fehlberg2(precision::Type{T} = Float64) where T <: AbstractFloat
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
    reconstructor = Fehlberg2

    return RungeKutta(name, butcher, iteration, reconstructor; ω)
end

"""
    Heun3(precision::Type{T} = Float64) where T <: AbstractFloat

Heun's third-order method.
"""
function Heun3(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Heun_3_2_1
    butcher = SMatrix{4, 6, precision, 24}(
        0, 0, 0, 0,
        1//3, 1//3, 0, 0,
        2//3, 0, 2//3, 0,
        1, 1//4, 0, 3//4,
        1, 0, 1//2, 1//2,
        1, 1, 0, 0
    ) |> transpose

    ω = SMatrix{3, 3, precision, 9}(
        1, -9//4, 3//2,
        0, 3, -3,
        0, -3//4, 3//2
    ) |> transpose

    iteration = Explicit()
    reconstructor = Heun3

    return RungeKutta(name, butcher, iteration, reconstructor; ω)
end

"""
    Ralston3(precision::Type{T} = Float64) where T <: AbstractFloat

Ralston's third-order method.

https://www.ams.org/journals/mcom/1962-16-080/S0025-5718-1962-0150954-0/S0025-5718-1962-0150954-0.pdf
"""
function Ralston3(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Ralston_3_2_1
    butcher = SMatrix{4, 6, precision, 24}(
        0, 0, 0, 0,
        1//2, 1//2, 0, 0,
        3//4, 0, 3//4, 0,
        1, 2//9, 1//3, 4//9,
        1, 1//3, 0, 2//3,
        1, 1, 0, 0
    ) |> transpose

    ω = SMatrix{3, 3, precision, 9}(
        1, -5//3, 8//9,
        0, 3, -8//3,
        0, -4//3, 16//9
    ) |> transpose

    iteration = Explicit()
    reconstructor = Ralston3

    return RungeKutta(name, butcher, iteration, reconstructor; ω)
end

"""
    Kutta3(precision::Type{T} = Float64) where T <: AbstractFloat

Kutta's third-order method.
"""
function Kutta3(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Kutta_3_2_1
    butcher = SMatrix{4, 6, precision, 24}(
        0, 0, 0, 0,
        1//2, 1//2, 0, 0,
        1, -1, 2, 0,
        1, 1//6, 2//3, 1//6,
        1, 1//4, 1//2, 1//4,
        1, 1, 0, 0
    ) |> transpose

    ω = SMatrix{3, 3, precision, 9}(
        1, -3//2, 2//3,
        0, 2, -4//3,
        0, -1//2, 2//3
    ) |> transpose

    iteration = Explicit()
    reconstructor = Kutta3

    return RungeKutta(name, butcher, iteration, reconstructor; ω)
end

"""
    BogackiShampine3(precision::Type{T} = Float64) where T <: AbstractFloat

Bogacki and Shampine's third-order method.

https://www.sciencedirect.com/science/article/pii/0893965989900797
"""
function BogackiShampine3(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Bogacki_Shampine_3_2_1
    butcher = SMatrix{5, 7, precision, 35}(
        0, 0, 0, 0, 0,
        1//2, 1//2, 0, 0, 0,
        3//4, 0, 3//4, 0, 0,
        1, 2//9, 1//3, 4//9, 0,
        1, 2//9, 1//3, 4//9, 0,
        1, 7//24, 1//4, 1//3, 1//8,
        1, 1, 0, 0, 0
    ) |> transpose

    # note: 3rd order C1 interpolant makes use of FSAL property
    ω = SMatrix{3, 4, precision, 12}(
        1, -4//3, 5//9,
        0, 1, -2//3,
        0, 4//3, -8//9,
        0, -1, 1
    ) |> transpose

    iteration = Explicit()
    reconstructor = BogackiShampine3

    return RungeKutta(name, butcher, iteration, reconstructor; ω)
end

"""
    ShuOsher3(precision::Type{T} = Float64) where T <: AbstractFloat

Shu and Osher's third-order SSP method.

https://www.sciencedirect.com/science/article/pii/0021999188901775
"""
function ShuOsher3(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Shu_Osher_3_2_1
    butcher = SMatrix{4, 6, precision, 24}(
        0, 0, 0, 0,
        1, 1, 0, 0,
        1//2, 1//4, 1//4, 0,
        1, 1//6, 1//6, 2//3,
        1, 4//9, 4//9, 1//9,
        1, 1, 0, 0
    ) |> transpose

    # note: only 2nd-order C0 interpolant preserves SSP
    ω = SMatrix{2, 3, precision, 6}(
        1, -5//6,
        0, 1//6,
        0, 2//3,
    ) |> transpose

    iteration = Explicit()
    reconstructor = ShuOsher3

    return RungeKutta(name, butcher, iteration, reconstructor; ω)
end

"""
    SpiteriRuuth3(precision::Type{T} = Float64) where T <: AbstractFloat

Spiteri and Ruuth's third-order SSP method.

https://epubs.siam.org/doi/10.1137/S0036142902419284
"""
function SpiteriRuuth3(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Spiteri_Ruuth_3_2_1
    butcher = SMatrix{5, 7, precision, 35}(
        0, 0, 0, 0, 0,
        1//2, 1//2, 0, 0, 0,
        1, 1//2, 1//2, 0, 0,
        1//2, 1//6, 1//6, 1//6, 0,
        1, 1//6, 1//6, 1//6, 1//2,
        1, 1//3, 1//6, 1//3, 1//6,
        1, 1, 0, 0, 0,
    ) |> transpose

    # note: only 2nd-order C0 interpolant preserves SSP
    ω = SMatrix{2, 4, precision, 8}(
        1, -5//6,
        0, 1//6,
        0, 1//6,
        0, 1//2,
    ) |> transpose

    iteration = Explicit()
    reconstructor = SpiteriRuuth3

    return RungeKutta(name, butcher, iteration, reconstructor; ω)
end
