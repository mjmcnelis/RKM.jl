
function reconstruct_method(method::RungeKutta,
                            precision::Type{T}) where T <: AbstractFloat

    return method.reconstructor(precision)
end

function reconstruct_method(method::LinearMultistep,
                            precision::Type{T}) where T <: AbstractFloat
    order = method.order
    start_method = method.start_method
    start_method = reconstruct_method(start_method, precision)

    # TODO: still not type-stable
    return method.reconstructor(; order, precision, start_method)
end

"""
    reconstruct_method(method::ODEMethod, precision::Type{T}) where T <: AbstractFloat

Reconstructs the ODE method (i.e. tables) with the prescribed numerical precision.

Required parameters: `method`, `precision`
"""
function reconstruct_method(method::ODEMethod, precision::Type{T}) where T <: AbstractFloat
    @warn "reconstruct_method(::$method ::$precision) not implemented!"
end
