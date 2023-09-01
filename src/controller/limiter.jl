
abstract type LimiterMethod end

# TODO: make piecewise and smooth limiter methods
# eventually these would also impose bounds by dt_min, dt_max
# if I re-evaluate rescale based on dt_min <= dt <= dt_max,
# would floating precision errors be a problem?

struct PiecewiseLimiter <: LimiterMethod
    """Safety factor to scale down estimate for predicted time step"""
    safety::Float64
    """Lower bound on the time step's rate of change"""
    low::Float64
    """Upper bound on the time step's rate of change"""
    high::Float64
end

function PiecewiseLimiter(; safety = 0.8, low = 0.2, high = 5.0)
    check_limiter_parameters(; safety, low, high)

    return PiecewiseLimiter(safety, low, high)
end

function check_limiter_parameters(; safety, low, high)
    @assert 0.0 <= low < 1.0 "low = $low is out of bounds [0, 1)"
    @assert 1.0 < high <= Inf "high = $high is out of bounds (1, Inf]"
    @assert 0.0 < safety < 1.0 "safety = $safety is out of bounds (0, 1)"
    @assert safety*high > 1.0 "safety*high = $(safety*high) is not greater than 1"
    return nothing
end

function limit_time_step(limiter::PiecewiseLimiter, rescale)
    @unpack safety, low, high = limiter
    return min(high, max(low, safety*rescale))
end