
function reconstruct_method(method::RungeKutta,
                            precision::Type{T}) where T <: AbstractFloat
    @unpack reconstructor = method
    return reconstructor(precision)
end

function reconstruct_method(method::LinearMultistep,
                            precision::Type{T}) where T <: AbstractFloat
    @unpack order, start_method = method

    start_method = reconstruct_method(start_method, precision)

    name = replace(String(method.name), "_" => "")
    name = filter(!isdigit, collect(name)) |> String |> Symbol
    return getfield(RKM, name)(; order, precision, start_method)
end

"""
    reconstruct_method(method::ODEMethod, precision::Type{T}) where T <: AbstractFloat

Reconstructs the ODE method (i.e. tables) with the prescribed numerical precision.

Required parameters: `method`, `precision`
"""
function reconstruct_method(method::ODEMethod, precision::Type{T}) where T <: AbstractFloat
    @warn "reconstruct_method(::$method ::$precision) not implemented!"
end
