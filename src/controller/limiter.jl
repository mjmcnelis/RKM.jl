
abstract type LimiterMethod end

function check_limiter_parameters(; safety, low, high, dt_min, dt_max)
    @assert 0.0 < low < 1.0 "low = $low is out of bounds (0, 1)"
    @assert 1.0 < high <= Inf "high = $high is out of bounds (1, Inf]"
    @assert 0.0 < safety < 1.0 "safety = $safety is out of bounds (0, 1)"
    @assert safety*high > 1.0 "safety*high = $(safety*high) is not greater than 1"
    @assert dt_min > 0 && dt_max > 0 "dt_min, dt_max = ($dt_min, $dt_max) are not positive"
    @assert dt_max > dt_min "dt_max = $dt_max is less than dt_min = $dt_min"
    return nothing
end

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
    """Minimum time step"""
    dt_min::Float64
    """Maximum time step"""
    dt_max::Float64
end

function PiecewiseLimiter(; safety = 0.8, low = 0.2, high = 5.0,
                            dt_min = eps(1.0), dt_max = Inf)
    check_limiter_parameters(; safety, low, high, dt_min, dt_max)

    return PiecewiseLimiter(safety, low, high, dt_min, dt_max)
end

struct SmoothLimiter <: LimiterMethod
    a::Float64
    b::Float64
    c::Float64
    d::Float64
    """Safety factor to scale down estimate for predicted time step"""
    safety::Float64
    """Lower bound on the time step's rate of change"""
    low::Float64
    """Upper bound on the time step's rate of change"""
    high::Float64
    """Minimum time step"""
    dt_min::Float64
    """Maximum time step"""
    dt_max::Float64
end

function SmoothLimiter(; safety = 0.8, low = 0.2, high = 5.0,
                         dt_min = eps(1.0), dt_max = Inf)
    check_limiter_parameters(; safety, low, high, dt_min, dt_max)

    # TODO: can fall back to 1 + (1-low)*tanh(log(x)/(1-low)) for large values of low
    if low > 0.6275
        @warn "Smooth limiter function is not monotonic for low = $low > 0.6275"
    end

    d = (low - 0.5) / (low - 0.5 - exp(-1.0/low)*(low + 0.5))
    a = low * (1.0 - d)
    b = d
    c = 0.5 * (1.0 - d*(1.0 - exp(-1.0/low)))

    return SmoothLimiter(a, b, c, d, safety, low, high, dt_min, dt_max)
end

function limit_time_step(limiter::PiecewiseLimiter, rescale)
    @unpack safety, low, high = limiter
    return min(high, max(low, safety*rescale))
end

function limit_time_step(limiter::SmoothLimiter, rescale)
    @unpack a, b, c, d, safety, low, high = limiter

    x = safety*rescale

    if x < 0.0
        @warn "rescale = $rescale is negative, limiting rescale = low"
        return low
    elseif 0.0 <= x <= 1.0
        return a + b*x + c*x^2 + d*low*exp(-x/low)

        # alternative for low > 0.6275
        # return 1.0 + (1.0 - low)*tanh(log(x + eps(x))/(1.0 - low))
    else
        return 1.0 + (high - 1.0)*tanh((x - 1.0)/(high - 1.0))
    end
end