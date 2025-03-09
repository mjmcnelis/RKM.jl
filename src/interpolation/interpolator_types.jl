
abstract type Interpolator end
abstract type DenseInterpolator <: Interpolator end

struct NoInterpolation <: Interpolator end
struct CubicHermite <: DenseInterpolator end
struct ContinuousFormula <: DenseInterpolator end

function interpolate_solution(interpolator::NoInterpolation, sol::Solution,
                              args...; kwargs...)
    @warn "No dense output for $(typeof(interpolator)), getting original solution..."
    return get_solution(sol)
end