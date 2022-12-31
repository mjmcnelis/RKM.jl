
abstract type AdaptiveStepSize end

struct Fixed <: AdaptiveStepSize end

struct FiniteDiff <: AdaptiveStepSize end
struct Doubling <: AdaptiveStepSize
    epsilon::Float64
    low::Float64
    high::Float64
    safety::Float64
    p_norm::Float64
    dt_min::Float64
    dt_max::Float64
    max_attempts::Int64
end

function Doubling(; epsilon = 1e-6, low = 0.2, high = 5.0, safety = 0.9, p_norm = 2,
                    dt_min = eps(1.0), dt_max = Inf, max_attempts = 10)

    check_adaptive_parameters(epsilon, low, high, safety, p_norm,
                              dt_min, dt_max, max_attempts)

    Doubling(epsilon, low, high, safety, p_norm, dt_min, dt_max, max_attempts)
end

struct Embedded <: AdaptiveStepSize
    epsilon::Float64
    low::Float64
    high::Float64
    safety::Float64
    p_norm::Float64
    dt_min::Float64
    dt_max::Float64
    max_attempts::Int
end

function Embedded(; epsilon = 1e-6, low = 0.2, high = 5.0, safety = 0.9, p_norm = 2,
                    dt_min = eps(1.0), dt_max = Inf, max_attempts = 10)

    check_adaptive_parameters(epsilon, low, high, safety, p_norm,
                              dt_min, dt_max, max_attempts)

    Embedded(epsilon, low, high, safety, p_norm, dt_min, dt_max, max_attempts)
end

function check_adaptive_parameters(epsilon, low, high, safety, p_norm,
                                   dt_min, dt_max, max_attempts)

    @assert epsilon > 0.0 "epsilon = $epsilon is not positive"
    @assert 0.0 <= low < 1.0 "low = $low is out of bounds [0, 1)"
    @assert 1.0 < high <= Inf "high = $high is out of bounds (1, Inf]"
    @assert 0.0 < safety < 1.0 "safety = $safety is out of bounds (0, 1)"
    @assert safety*high > 1.0 "safety*high = $(safety*high) is not greater than 1"
    @assert p_norm >= 1 "p_norm = $p_norm is not valid"
    @assert dt_min > 0 && dt_max > 0 "dt_min, dt_max = ($dt_min, $dt_max) are not positive"
    @assert dt_max > dt_min "dt_max = $dt_max is less than dt_min = $dt_min"
    @assert max_attempts > 0 "max_attempts = $max_attempts is not greater than 0"
end
