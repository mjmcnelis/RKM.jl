
abstract type AdaptiveStepSize end

struct Fixed <: AdaptiveStepSize end

struct CentralDiff <: AdaptiveStepSize
    """Relative and incremental error tolerance"""
    epsilon::Float64
    """Integer used to compute L--norms"""
    p_norm::Float64
end

function CentralDiff(; epsilon = 1e-6, p_norm = 2)

    check_adaptive_parameters_1(; epsilon, p_norm)

    return CentralDiff(epsilon, p_norm)
end

"""
$(TYPEDEF)

Step doubling adaptive time step algorithm

# Fields
$(TYPEDFIELDS)
"""
struct Doubling <: AdaptiveStepSize
    """Relative error tolerance"""
    epsilon::Float64
    """Absolute error tolerance"""
    alpha::Float64
    """Incremental error tolerance"""
    delta::Float64
    """Integer used to compute L--norms"""
    p_norm::Float64
    """Maximum number of attempts to compute time step per update"""
    max_attempts::Int64
    """Total number of attempts in evolution loop"""
    total_attempts::MVector{1,Int64}
end

function Doubling(; epsilon = 1e-6, alpha = 1e-6, delta = 1e-6,
                    p_norm = 2, max_attempts = 10)

    check_adaptive_parameters_1(; epsilon, alpha, delta, p_norm)
    check_adaptive_parameters_2(; max_attempts)
    total_attempts = MVector{1,Int64}(0)

    return Doubling(epsilon, alpha, delta, p_norm, max_attempts, total_attempts)
end

struct Embedded <: AdaptiveStepSize
    """Relative error tolerance"""
    epsilon::Float64
    """Absolute error tolerance"""
    alpha::Float64
    """Incremental error tolerance"""
    delta::Float64
    """Integer used to compute L--norms"""
    p_norm::Float64
    """Maximum number of attempts to compute time step per update"""
    max_attempts::Int64
    """Total number of attempts in evolution loop"""
    total_attempts::MVector{1,Int64}
end

function Embedded(; epsilon = 1e-6, alpha = 1e-6, delta = 1e-6,
                    p_norm = 2, max_attempts = 10)

    check_adaptive_parameters_1(; epsilon, alpha, delta, p_norm)
    check_adaptive_parameters_2(; max_attempts)
    total_attempts = MVector{1,Int64}(0)

    return Embedded(epsilon, alpha, delta, p_norm, max_attempts, total_attempts)
end

function check_adaptive_parameters_1(; epsilon, alpha, delta, p_norm)
    @assert epsilon >= 0.0 "epsilon = $epsilon cannot be negative"
    @assert alpha >= 0.0 "alpha = $alpha cannot be negative"
    @assert delta >= 0.0 "delta = $delta cannot be negative"
    msg = "one of the tolerance parameters (epsilon, alpha, delta) must be positive"
    @assert !all(x -> x == 0.0, [epsilon, alpha, delta]) msg
    @assert p_norm >= 1 "p_norm = $p_norm is not valid"
    return nothing
end

function check_adaptive_parameters_2(; max_attempts)
    @assert max_attempts > 0 "max_attempts = $max_attempts is not greater than 0"
    return nothing
end

function reset_attempts!(adaptive::AdaptiveStepSize)
    if hasproperty(adaptive, :total_attempts)
        adaptive.total_attempts[1] = 0
    end
    return nothing
end

function compute_step_rejection_rate(adaptive::AdaptiveStepSize, timer::TimeLimit)
    if hasproperty(adaptive, :total_attempts)
        @unpack total_steps = timer
        @unpack total_attempts = adaptive
        return 100.0 * (1.0 - total_steps[1]/total_attempts[1])
    end
    return 0.0
end