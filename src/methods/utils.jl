
function reconstruct_method(method::RungeKutta, precision::Type{T}) where T <: AbstractFloat
    name = replace(String(method.name), "_" => "") |> Symbol
    return getfield(RKM, name)(; precision)
end

function reconstruct_method(method::Adams, precision::Type{T}) where T <: AbstractFloat
    @unpack order = method
    name = replace(String(method.name), "_" => "")
    name = filter(!isdigit, collect(name)) |> String |> Symbol
    return getfield(RKM, name)(; order, precision)
end

"""
    reconstruct_method(method::ODEMethod, precision::Type{T}) where T <: AbstractFloat

Reconstructs the ODE method (i.e. tables) with the prescribed numerical precision.

Required parameters: `method`, `precision`
"""
function reconstruct_method(method::ODEMethod, precision::Type{T}) where T <: AbstractFloat
    @warn "reconstruct_method(::$method ::$precision) not implemented!"
end
