"""
    reconstruct_method(method::ODEMethod, precision::Type{T}) where T <: AbstractFloat

Reconstructs the ODE method (i.e. tables) with the prescribed numerical precision.

Required parameters: `method`, `precision`
"""
function reconstruct_method(method::ODEMethod, precision::Type{T}) where T <: AbstractFloat
    name = replace(String(method.name), "_" => "") |> Symbol
    return getfield(RKM, name)(; precision)
end