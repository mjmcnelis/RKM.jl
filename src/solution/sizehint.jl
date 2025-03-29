"""
    sizehint_solution!(adaptive::Fixed, interpolator::Interpolator, sol::Solution,
                       t0::T, tf::T, dt::T, sensitivity::SensitivityMethod,
                       save_time_derivative::Bool,
                       stages::Int64) where T <: AbstractFloat

Applies `sizehint!` to the time `t` and state variables `y` in the solution `sol`.
The time derivative `f`, intermediate stages `dy` and sensitivity coefficients `S`
are also size-hinted if they are being outputted.

Required parameters: `adaptive`, `interpolator`, `sol`, `t0`, `tf`, `dt`,
                     `sensitivity`, `save_time_derivative`, `stages`
"""
function sizehint_solution!(adaptive::Fixed, interpolator::Interpolator, sol::Solution,
                            t0::T, tf::T, dt::T, sensitivity::SensitivityMethod,
                            save_time_derivative::Bool,
                            stages::Int64) where T <: AbstractFloat
    ny = sol.dimensions[1]
    np = sol.coefficients[1]
    # note: add 2 instead of 1 b/c of round-off errors
    nt = 2 + round(Int64, Float64((tf - t0)/dt))

    sizehint!(sol.t, nt)
    sizehint!(sol.y, ny*nt)
    if save_time_derivative || interpolator isa CubicHermite
        sizehint!(sol.f, ny*nt)
    end
    if interpolator isa ContinuousFormula
        # TODO: change back to nt if decide to output dy at final time
        sizehint!(sol.dy, ny*(nt-1)*stages)
    end
    if !(sensitivity isa NoSensitivity)
        sizehint!(sol.S, ny*np*nt)
    end
    return nothing
end

function sizehint_solution!(adaptive::AdaptiveStepSize, args...)
    return nothing
end
