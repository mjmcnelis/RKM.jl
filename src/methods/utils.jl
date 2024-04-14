
function reconstruct_method(method::RungeKutta,
                            precision::Type{T}) where T <: AbstractFloat

    @unpack name, c, A_T, b, b_hat, stages, order,
    iteration, fsal, explicit_stage, code_name = method

    # ugh this won't fix precision (unless default is BigFloat)

    # convert to type precision (maybe need to do this w/ other reconstruct functions)
    c2 = precision.(c)
    A_T2 = precision.(A_T)
    b2 = precision.(b)
    b_hat2 = precision.(b_hat)
    order2 = precision.(order)

    return RungeKutta(name, c2, A_T2, b2, b_hat2, stages, order2,
                      iteration, fsal, explicit_stage, code_name)
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
