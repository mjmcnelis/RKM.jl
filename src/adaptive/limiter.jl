
abstract type LimiterMethod end

@kwdef struct PiecewiseLimiter <: LimiterMethod
    """Safety factor to scale down estimate for predicted time step"""
    safety::Float64 = 0.8
    """Lower bound on the time step's rate of change"""
    low::Float64 = 0.2
    """Upper bound on the time step's rate of change"""
    high::Float64 = 5.0
    """Minimum time step"""
    dt_min::Float64 = eps(1.0)
    """Maximum time step"""
    dt_max::Float64 = Inf

    function PiecewiseLimiter(safety, low, high, dt_min, dt_max)
        st = new(safety, low, high, dt_min, dt_max)
        check_limiter(st)
        return st
    end
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

    function SmoothLimiter(a, b, c, d, safety, low, high, dt_min, dt_max)
        st = new(a, b, c, d, safety, low, high, dt_min, dt_max)
        check_limiter(st)
        return st
    end
end

function SmoothLimiter(; safety = 0.8, low = 0.2, high = 5.0,
                         dt_min = eps(1.0), dt_max = Inf)
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

function check_limiter(limiter::LM) where LM <: LimiterMethod
    safety = limiter.safety
    low = limiter.low
    high = limiter.high
    dt_min = limiter.dt_min
    dt_max = limiter.dt_max

    @assert 0.0 < low < 1.0 "low = $low is out of bounds (0, 1)"
    @assert 1.0 < high <= Inf "high = $high is out of bounds (1, Inf]"
    @assert 0.0 < safety < 1.0 "safety = $safety is out of bounds (0, 1)"
    # note: need following assert b/c I set rescale = high when error = 0
    @assert safety*high > 1.0 "safety*high = $(safety*high) is not greater than 1"
    @assert dt_min > 0 && dt_max > 0 "dt_min, dt_max = ($dt_min, $dt_max) are not positive"
    @assert dt_max > dt_min "dt_max = $dt_max is less than dt_min = $dt_min"

    return nothing
end

function limit_time_step(limiter::PiecewiseLimiter, rescale::T) where T <: AbstractFloat
    safety = limiter.safety
    low = limiter.low
    high = limiter.high

    return min(high, max(low, safety*rescale))
end

function limit_time_step(limiter::SmoothLimiter, rescale::T) where T <: AbstractFloat
    a = limiter.a
    b = limiter.b
    c = limiter.c
    d = limiter.d
    safety = limiter.safety
    low = limiter.low
    high = limiter.high

    x = T(safety)*rescale

    if x < 0.0
        return T(low)
    elseif 0.0 <= x <= 1.0
        return a + b*x + c*x^2 + d*low*exp(-x/low)

        # alternative for low > 0.6275
        # return 1.0 + (1.0 - low)*tanh(log(x + eps(x))/(1.0 - low))
    else
        return 1.0 + (high - 1.0)*tanh((x - 1.0)/(high - 1.0))
    end
end