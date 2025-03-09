"""
    sizehint_solution!(adaptive::Fixed, interpolator::Interpolator,
                       sol::Solution, t0::T, tf::T, dt::T,
                       sensitivity_method::SensitivityMethod,
                       save_time_derivative::Bool) where T <: AbstractFloat

Applies `sizehint!` to the state vector `y` and time `t` in the solution `sol`.
The ODE function vector `f` and sensitivity coefficients `S` are also size-hinted
if they are being outputted.

Required parameters: `adaptive`, `interpolator`, `sol`, `t0`, `tf`,
                     `dt`, `sensitivity_method`, `save_time_derivative`
"""
function sizehint_solution!(adaptive::Fixed, interpolator::Interpolator,
                            sol::Solution, t0::T, tf::T, dt::T,
                            sensitivity_method::SensitivityMethod,
                            save_time_derivative::Bool) where T <: AbstractFloat
    ny = sol.dimensions[1]
    np = sol.coefficients[1]
    nt = 2 + round(Int64, Float64((tf - t0)/dt))

    sizehint!(sol.t, nt)
    sizehint!(sol.y, ny*nt)
    if save_time_derivative || interpolator isa DenseInterpolator
        sizehint!(sol.f, ny*nt)
    end
    if !(sensitivity_method isa NoSensitivity)
        sizehint!(sol.S, ny*np*nt)
    end
    return nothing
end

function sizehint_solution!(adaptive::AdaptiveStepSize, args...)
    return nothing
end
