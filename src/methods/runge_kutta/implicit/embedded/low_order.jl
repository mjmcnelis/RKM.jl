"""
    CrankNicolson21(; precision::Type{T} = Float64) where T <: AbstractFloat

Crank and Nicolson's second-order method (A-stable) with embedded Euler pair.
"""
function CrankNicolson21(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [0 0 0
               1 1//2 1//2
               1 1//2 1//2
               1 1 0]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Crank_Nicolson_2_1, butcher)
end

function LobattoIIIB21(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [1//2 1//2 0
               1//2 1//2 0
               1 1//2 1//2
               1 1 0]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Lobatto_IIIB_2_1, butcher)
end

function LobattoIIIC21(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [0 1//2 -1//2
               1 1//2 1//2
               1 1//2 1//2
               1 1 0]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Lobatto_IIIC_2_1, butcher)
end
