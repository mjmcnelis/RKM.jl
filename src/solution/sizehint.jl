
"""
    _sizehint_solution!(sol::Solution, t0::T, tf::T, dt::T1,
         sensitivity_method::SensitivityMethod) where {T <: AbstractFloat,
                                                       T1 <: AbstractFloat}

Applies `sizehint!` to the vector fields `y` and `t` in the solution `sol`.

Required parameters: `sol`, `t0`, `tf`, `dt0`, `sensitivity_method`
"""
function _sizehint_solution!(sol::Solution, t0::T, tf::T, dt::T1,
             sensitivity_method::SensitivityMethod) where {T <: AbstractFloat,
                                                           T1 <: AbstractFloat}
    dimensions = sol.dimensions[1]
    coefficients = sol.coefficients[1]

    steps = round(Int64, Float64((tf - t0)/dt))
    # TODO: can I use steps + 1 if dense output?
    sizehint!(sol.y, dimensions*(steps + 2))
    sizehint!(sol.t, steps + 2)
    if !(sensitivity_method isa NoSensitivity)
        sizehint!(sol.S, dimensions*coefficients*(steps + 2))
    end
    return nothing
end

function sizehint_solution!(::AdaptiveStepSize, ::NoInterpolator, args...)
    return nothing
end

function sizehint_solution!(::Fixed, ::NoInterpolator, args...)
    # note: may take (steps+1) steps before solver finishes b/c of float precision errors
    return _sizehint_solution!(args...)
end

function sizehint_solution!(::AdaptiveStepSize, interpolator::DenseInterpolator,
                            sol::Solution, t0::T, tf::T, dt0::T,
                            sensitivity_method::SensitivityMethod) where T <: AbstractFloat
    @unpack dt_save = interpolator
    return _sizehint_solution!(sol, t0, tf, dt_save, sensitivity_method)
end
