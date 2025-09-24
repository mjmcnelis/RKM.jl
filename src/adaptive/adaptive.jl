
abstract type AdaptiveTimeStep end

struct Fixed <: AdaptiveTimeStep end

@kwdef struct CentralDiff{LM} <: AdaptiveTimeStep where LM <: LimiterMethod
    """Relative error tolerance"""
    epsilon::Float64 = 1e-6
    """Absolute error tolerance"""
    alpha::Float64 = 1e-6
    """Incremental error tolerance"""
    delta::Float64 = 1e-6
    """Integer used to compute L--norms"""
    p_norm::Float64 = 2.0
    """Skip rescaling tolerance parameters if benchmark against OrdinaryDiffEq"""
    benchmark_diffeq::Bool = false
    """Limiter method for time step controller"""
    limiter::LM = PiecewiseLimiter()

    function CentralDiff(epsilon, alpha, delta, p_norm, benchmark_diffeq, limiter)
        st = new{typeof(limiter)}(epsilon, alpha, delta, p_norm, benchmark_diffeq, limiter)
        check_adaptive(st)
        return st
    end
end

"""
$(TYPEDEF)

Step doubling adaptive time step algorithm

# Fields
$(TYPEDFIELDS)
"""
@kwdef struct Doubling{PCB, LM} <: AdaptiveTimeStep where {PCB <: PIDControlBeta,
                                                           LM <: LimiterMethod}
    """Relative error tolerance"""
    epsilon::Float64 = 1e-6
    """Absolute error tolerance"""
    alpha::Float64 = 1e-6
    """Incremental error tolerance"""
    delta::Float64 = 1e-6
    """Integer used to compute L--norms"""
    p_norm::Float64 = 2.0
    """Maximum number of attempts to compute time step per update"""
    max_attempts::Int64 = 10
    """Total number of attempts in evolution loop"""
    total_attempts::MVector{1,Int64} = MVector{1,Int64}(0)
    """Skip rescaling tolerance parameters if benchmark against OrdinaryDiffEq"""
    benchmark_diffeq::Bool = false
    """PID controller gain parameters"""
    pid::PCB = PIControl()
    """Limiter method for time step controller"""
    limiter::LM = PiecewiseLimiter()
    """Whether the time step controller has been initialized"""
    initialized_controller::MVector{1,Bool} = MVector{1, Bool}(false)

    function Doubling(epsilon, alpha, delta, p_norm, max_attempts, total_attempts,
                      benchmark_diffeq, pid, limiter, initialized_controller)
        st = new{typeof(pid), typeof(limiter)}(epsilon, alpha, delta, p_norm,
                                               max_attempts, total_attempts,
                                               benchmark_diffeq, pid, limiter,
                                               initialized_controller)
        check_adaptive(st)
        return st
    end
end

@kwdef struct Embedded{PCB, LM} <: AdaptiveTimeStep where {PCB <: PIDControlBeta,
                                                           LM <: LimiterMethod}
    """Relative error tolerance"""
    epsilon::Float64 = 1e-6
    """Absolute error tolerance"""
    alpha::Float64 = 1e-6
    """Incremental error tolerance"""
    delta::Float64 = 1e-6
    """Integer used to compute L--norms"""
    p_norm::Float64 = 2.0
    """Maximum number of attempts to compute time step per update"""
    max_attempts::Int64 = 10
    """Total number of attempts in evolution loop"""
    total_attempts::MVector{1,Int64} = MVector{1,Int64}(0)
    """Skip rescaling tolerance parameters if benchmark against OrdinaryDiffEq"""
    benchmark_diffeq::Bool = false
    """PID controller gain parameters"""
    pid::PCB = PIControl()
    """Limiter method for time step controller"""
    limiter::LM = PiecewiseLimiter()
    """Whether the time step controller has been initialized"""
    initialized_controller::MVector{1,Bool} = MVector{1, Bool}(false)

    function Embedded(epsilon, alpha, delta, p_norm, max_attempts, total_attempts,
                      benchmark_diffeq, pid, limiter, initialized_controller)
        st = new{typeof(pid), typeof(limiter)}(epsilon, alpha, delta, p_norm,
                                               max_attempts, total_attempts,
                                               benchmark_diffeq, pid, limiter,
                                               initialized_controller)
        check_adaptive(st)
        return st
    end
end

function check_adaptive(adaptive::Fixed)
    return nothing
end

function check_adaptive(adaptive::AdaptiveTimeStep)

    epsilon = adaptive.epsilon
    alpha = adaptive.alpha
    delta = adaptive.delta
    p_norm = adaptive.p_norm

    @assert epsilon >= 0.0 "epsilon = $epsilon cannot be negative"
    @assert alpha >= 0.0 "alpha = $alpha cannot be negative"
    @assert delta >= 0.0 "delta = $delta cannot be negative"
    msg = "one of the tolerance parameters (epsilon, alpha, delta) must be positive"
    @assert !all(x -> x == 0.0, [epsilon, alpha, delta]) msg
    @assert p_norm >= 1 "p_norm = $p_norm is not valid"

    if hasproperty(adaptive, :max_attempts)
        max_attempts = adaptive.max_attempts
        @assert max_attempts > 0 "max_attempts = $max_attempts is not greater than 0"
    end

    return nothing
end

function get_local_order(adaptive::CentralDiff,
                         order::SVector{P,T}) where {P, T <: AbstractFloat}
    return 2.0
end

function get_local_order(adaptive::Doubling,
                         order::SVector{P,T}) where {P, T <: AbstractFloat}
    return order[1] + 1.0
end

function get_local_order(adaptive::Embedded,
                         order::SVector{P,T}) where {P, T <: AbstractFloat}
    # TODO: should be more dynamic, depend on largest error pair
    #       for now just assume first embedded pair
    return order[2] + 1.0
end

function rescale_tolerance(adaptive::CentralDiff,
                           order::SVector{P,T}) where {P, T <: AbstractFloat}
    return 1.0 / order[1]
end

function rescale_tolerance(adaptive::Doubling,
                           order::SVector{P,T}) where {P, T <: AbstractFloat}
    return order[1] / (1.0 + order[1])
end

function rescale_tolerance(adaptive::Embedded,
                           order::SVector{P,T}) where {P, T <: AbstractFloat}
    # TODO: should be more dynamic, depend on largest error pair
    #       for now just assume first embedded pair
    return order[2] / order[1]
end

function reconstruct_adaptive(adaptive::Fixed,
                              order::SVector{P,T}) where {P, T <: AbstractFloat}
    return adaptive
end

function reconstruct_adaptive(adaptive::Fixed, order::Int64)
    return adaptive
end

function reconstruct_adaptive(adaptive::AdaptiveTimeStep,
                              order::SVector{P,T}) where {P, T <: AbstractFloat}

    benchmark_diffeq = adaptive.benchmark_diffeq
    limiter = adaptive.limiter

    if !(adaptive isa CentralDiff)
        @set! adaptive.total_attempts = MVector{1,Int64}(0)
        @set! adaptive.initialized_controller = MVector{1,Bool}(false)
    end

    # rescale tolerance parameters
    if !benchmark_diffeq
        repower_high = rescale_tolerance(adaptive, order) |> Float64

        @set! adaptive.epsilon ^= repower_high
        @set! adaptive.alpha ^= repower_high
        @set! adaptive.delta ^= repower_high
        @set! limiter.high ^= repower_high
        @set! adaptive.limiter = limiter
    end

    if !(adaptive isa CentralDiff)
        local_order = get_local_order(adaptive, order) |> Float64
        pid = adaptive.pid

        @set! pid.beta1 /= local_order
        @set! pid.beta2 /= local_order
        @set! pid.beta3 /= local_order
        @set! adaptive.pid = pid
    end

    check_adaptive(adaptive)

    return adaptive
end

"""
    compute_step_rejection_rate(adaptive::AdaptiveTimeStep, timer::TimeLimit)

Computes the step rejection rate [%] if the time step is adaptive.

Required parameters: `adaptive`, `timer`
"""
function compute_step_rejection_rate(adaptive::AdaptiveTimeStep, timer::TimeLimit)

    if hasproperty(adaptive, :total_attempts)
        total_steps = timer.total_steps
        total_attempts = adaptive.total_attempts

        return 100.0 * (1.0 - total_steps[1]/total_attempts[1])
    else
        return 0.0
    end
end