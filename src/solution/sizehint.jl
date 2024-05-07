
"""
    _sizehint_solution!(sol::Solution, t0::T, tf::T, dt::T1,
                        dimensions::Int64) where {T <: AbstractFloat,
                                                  T1 <: AbstractFloat}

Applies `sizehint!` to the vector fields `y` and `t` in the solution `sol`.

Required parameters: `sol`, `t0`, `tf`, `dt0`, `dimensions`
"""
function _sizehint_solution!(sol::Solution, t0::T, tf::T, dt::T1,
                             dimensions::Int64) where {T <: AbstractFloat,
                                                       T1 <: AbstractFloat}
    steps = round(Int64, Float64((tf - t0)/dt))
    # TODO: can I use steps + 1 if dense output?
    sizehint!(sol.y, dimensions*(steps + 2))
    sizehint!(sol.t, steps + 2)
    return nothing
end

sizehint_solution!(::AdaptiveStepSize, ::NoInterpolator, args...) = nothing

function sizehint_solution!(::Fixed, ::NoInterpolator, args...)
    # note: may take (steps+2) steps before solver finishes b/c of float precision errors
    return _sizehint_solution!(args...)
end

function sizehint_solution!(::AdaptiveStepSize, interpolator::DenseInterpolator,
                            sol::Solution, t0::T, tf::T, dt0::T,
                            dimensions::Int64) where T <: AbstractFloat
    @unpack dt_save = interpolator
    return _sizehint_solution!(sol, t0, tf, dt_save, dimensions)
end
