
abstract type AdaptiveStepSize end

struct Fixed <: AdaptiveStepSize end

struct CentralDiff <: AdaptiveStepSize
    """Relative and incremental error tolerance"""
    epsilon::Float64
    """Lower bound on the time step's rate of change"""
    low::Float64
    """Upper bound on the time step's rate of change"""
    high::Float64
    """Integer used to compute L--norms"""
    p_norm::Float64
    """Minimum time step"""
    dt_min::Float64
    """Maximum time step"""
    dt_max::Float64
end

function CentralDiff(; epsilon = 1e-6, low = 0.2, high = 1.5, p_norm = 2,
                      dt_min = eps(1.0), dt_max = Inf)

    check_adaptive_parameters_1(; epsilon, low, high, p_norm, dt_min, dt_max)

    return CentralDiff(epsilon, low, high, p_norm, dt_min, dt_max)
end

"""
$(TYPEDEF)

Step doubling adaptive time step algorithm

# Fields
$(TYPEDFIELDS)
"""
struct Doubling <: AdaptiveStepSize
    """Relative and incremental error tolerance"""
    epsilon::Float64
    """Lower bound on the time step's rate of change"""
    low::Float64
    """Upper bound on the time step's rate of change"""
    high::Float64
    """Safety factor to scale down estimate for new time step"""
    safety::Float64
    """Integer used to compute L--norms"""
    p_norm::Float64
    """Minimum time step"""
    dt_min::Float64
    """Maximum time step"""
    dt_max::Float64
    """Maximum number of attempts to search for new time step"""
    max_attempts::Int64
end

function Doubling(; epsilon = 1e-6, low = 0.2, high = 5.0, safety = 0.8, p_norm = 2,
                    dt_min = eps(1.0), dt_max = Inf, max_attempts = 10)

    check_adaptive_parameters_1(; epsilon, low, high, p_norm, dt_min, dt_max)
    check_adaptive_parameters_2(; high, safety, max_attempts)

    return Doubling(epsilon, low, high, safety, p_norm, dt_min, dt_max, max_attempts)
end

struct Embedded <: AdaptiveStepSize
    """Relative and incremental error tolerance"""
    epsilon::Float64
    """Lower bound on the time step's rate of change"""
    low::Float64
    """Upper bound on the time step's rate of change"""
    high::Float64
    """Safety factor to scale down estimate for new time step"""
    safety::Float64
    """Integer used to compute L--norms"""
    p_norm::Float64
    """Minimum time step"""
    dt_min::Float64
    """Maximum time step"""
    dt_max::Float64
    """Maximum number of attempts to search for new time step"""
    max_attempts::Int64
end

function Embedded(; epsilon = 1e-6, low = 0.2, high = 5.0, safety = 0.8, p_norm = 2,
                    dt_min = eps(1.0), dt_max = Inf, max_attempts = 10)

    check_adaptive_parameters_1(; epsilon, low, high, p_norm, dt_min, dt_max)
    check_adaptive_parameters_2(; high, safety, max_attempts)

    return Embedded(epsilon, low, high, safety, p_norm, dt_min, dt_max, max_attempts)
end

function check_adaptive_parameters_1(; epsilon, low, high, p_norm, dt_min, dt_max)
    @assert epsilon > 0.0 "epsilon = $epsilon is not positive"
    @assert 0.0 <= low < 1.0 "low = $low is out of bounds [0, 1)"
    @assert 1.0 < high <= Inf "high = $high is out of bounds (1, Inf]"
    @assert p_norm >= 1 "p_norm = $p_norm is not valid"
    @assert dt_min > 0 && dt_max > 0 "dt_min, dt_max = ($dt_min, $dt_max) are not positive"
    @assert dt_max > dt_min "dt_max = $dt_max is less than dt_min = $dt_min"
    return nothing
end

function check_adaptive_parameters_2(; high, safety, max_attempts)
    @assert 0.0 < safety < 1.0 "safety = $safety is out of bounds (0, 1)"
    @assert safety*high > 1.0 "safety*high = $(safety*high) is not greater than 1"
    @assert max_attempts > 0 "max_attempts = $max_attempts is not greater than 0"
    return nothing
end
