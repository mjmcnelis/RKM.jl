
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
    """Total number of attempts in evolution loop""" # TODO: move to update cache or solver
    total_attempts::MVector{1,Int64} = MVector{1,Int64}(0)
    """Skip rescaling tolerance parameters if benchmark against OrdinaryDiffEq"""
    benchmark_diffeq::Bool = false
    """PID controller gain parameters"""
    pid::PCB = PIControl()
    """Limiter method for time step controller"""
    limiter::LM = PiecewiseLimiter()
    """Whether the time step controller has been initialized"""
    initialized_controller::MVector{1,Bool} = MVector{1, Bool}(false)
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
    """Total number of attempts in evolution loop""" # TODO: move to update cache or solver
    total_attempts::MVector{1,Int64} = MVector{1,Int64}(0)
    """Skip rescaling tolerance parameters if benchmark against OrdinaryDiffEq"""
    benchmark_diffeq::Bool = false
    """PID controller gain parameters"""
    pid::PCB = PIControl()
    """Limiter method for time step controller"""
    limiter::LM = PiecewiseLimiter()
    """Whether the time step controller has been initialized"""
    initialized_controller::MVector{1,Bool} = MVector{1, Bool}(false)
end

function check_adaptive_parameters(adaptive::Fixed)
    return nothing
end

function check_adaptive_parameters(adaptive::AdaptiveTimeStep)

    epsilon = adaptive.epsilon
    alpha = adaptive.alpha
    delta = adaptive.delta
    p_norm = adaptive.p_norm
    max_attempts = adaptive.max_attempts

    @assert epsilon >= 0.0 "epsilon = $epsilon cannot be negative"
    @assert alpha >= 0.0 "alpha = $alpha cannot be negative"
    @assert delta >= 0.0 "delta = $delta cannot be negative"
    msg = "one of the tolerance parameters (epsilon, alpha, delta) must be positive"
    @assert !all(x -> x == 0.0, [epsilon, alpha, delta]) msg
    @assert p_norm >= 1 "p_norm = $p_norm is not valid"
    @assert max_attempts > 0 "max_attempts = $max_attempts is not greater than 0"

    return nothing
end

# function get_local_order(::CentralDiff, args...)
#     return 2.0
# end

function get_local_order(::Doubling, order::SVector{P,T}) where {P, T <: AbstractFloat}
    return order[1] + 1.0
end

function get_local_order(::Embedded, order::SVector{P,T}) where {P, T <: AbstractFloat}
    return minimum(order) + 1.0
end

# function rescale_tolerance(::CentralDiff, order::SVector{P,T}) where {P,T <: AbstractFloat}
#     return 1.0 / order[1]
# end

function rescale_tolerance(::Doubling, order::SVector{P,T}) where {P, T <: AbstractFloat}
    return order[1] / (1.0 + order[1])
end

function rescale_tolerance(::Embedded, order::SVector{P,T}) where {P, T <: AbstractFloat}
    return minimum(order) / maximum(order)
end

function reconstruct_adaptive(adaptive::Fixed, method::ODEMethod)
    return adaptive
end

function reconstruct_adaptive(adaptive::CentralDiff, method::ODEMethod)
    return adaptive
end

function reconstruct_adaptive(adaptive::AdaptiveTimeStep, method::ODEMethod)

    benchmark_diffeq = adaptive.benchmark_diffeq
    pid = adaptive.pid
    limiter = adaptive.limiter

    order = method.order

    @set! adaptive.total_attempts = MVector{1,Int64}(0)
    @set! adaptive.initialized_controller = MVector{1,Bool}(false)

    # rescale tolerance parameters
    if !benchmark_diffeq
        repower_high = rescale_tolerance(adaptive, order) |> Float64

        @set! adaptive.epsilon ^= repower_high
        @set! adaptive.alpha ^= repower_high
        @set! adaptive.delta ^= repower_high
        @set! limiter.high ^= repower_high
        @set! adaptive.limiter = limiter

        # TODO: add message
        # reminder: reason I need this is b/c I default rescale = high when error = 0
        # @code_warntype red %88 = Base.AssertionError("limiter.safety * limiter.high > 1.0")::Any
        @assert limiter.safety * limiter.high > 1.0
    end

    # TODO: why don't I need |> Float64 here...
    local_order = get_local_order(adaptive, order)

    @set! pid.beta1 /= local_order
    @set! pid.beta2 /= local_order
    @set! pid.beta3 /= local_order
    @set! adaptive.pid = pid

    check_adaptive_parameters(adaptive)

    return adaptive
end

function compute_step_rejection_rate(adaptive::AdaptiveTimeStep, timer::TimeLimit)

    if hasproperty(adaptive, :total_attempts)
        total_steps = timer.total_steps
        total_attempts = adaptive.total_attempts

        return 100.0 * (1.0 - total_steps[1]/total_attempts[1])
    else
        return 0.0
    end
end