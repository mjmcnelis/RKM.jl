
abstract type AdaptiveStepSize end

struct Fixed <: AdaptiveStepSize end

struct CentralDiff{T} <: AdaptiveStepSize where T <: AbstractFloat
    """Relative error tolerance"""
    epsilon::T
    """Absolute error tolerance"""
    alpha::T
    """Incremental error tolerance"""
    delta::T
    """Integer used to compute L--norms"""
    p_norm::T
    """Skip rescaling tolerance parameters if benchmark against OrdinaryDiffEq"""
    benchmark_diffeq::Bool
end

function CentralDiff(; epsilon = 1e-6, alpha = 1e-6, delta = 1e-6,
                       p_norm = 2.0, benchmark_diffeq = false)

    check_adaptive_parameters_1(; epsilon, alpha, delta, p_norm)

    return CentralDiff(epsilon, alpha, delta, p_norm, benchmark_diffeq)
end

"""
$(TYPEDEF)

Step doubling adaptive time step algorithm

# Fields
$(TYPEDFIELDS)
"""
struct Doubling{T} <: AdaptiveStepSize where T <: AbstractFloat
    """Relative error tolerance"""
    epsilon::T
    """Absolute error tolerance"""
    alpha::T
    """Incremental error tolerance"""
    delta::T
    """Integer used to compute L--norms"""
    p_norm::T
    """Maximum number of attempts to compute time step per update"""
    max_attempts::Int64
    """Total number of attempts in evolution loop"""
    total_attempts::MVector{1,Int64}
    """Skip rescaling tolerance parameters if benchmark against OrdinaryDiffEq"""
    benchmark_diffeq::Bool
end

function Doubling(; epsilon = 1e-6, alpha = 1e-6, delta = 1e-6,p_norm = 2.0,
                    max_attempts = 10, benchmark_diffeq = false)

    check_adaptive_parameters_1(; epsilon, alpha, delta, p_norm)
    check_adaptive_parameters_2(; max_attempts)
    total_attempts = MVector{1,Int64}(0)

    return Doubling(epsilon, alpha, delta, p_norm, max_attempts,
                    total_attempts, benchmark_diffeq)
end

struct Embedded{T} <: AdaptiveStepSize where T <: AbstractFloat
    """Relative error tolerance"""
    epsilon::T
    """Absolute error tolerance"""
    alpha::T
    """Incremental error tolerance"""
    delta::T
    """Integer used to compute L--norms"""
    p_norm::T
    """Maximum number of attempts to compute time step per update"""
    max_attempts::Int64
    """Total number of attempts in evolution loop"""
    total_attempts::MVector{1,Int64}
    """Skip rescaling tolerance parameters if benchmark against OrdinaryDiffEq"""
    benchmark_diffeq::Bool
end

function Embedded(; epsilon = 1e-6, alpha = 1e-6, delta = 1e-6, p_norm = 2.0,
                    max_attempts = 10, benchmark_diffeq = false)

    check_adaptive_parameters_1(; epsilon, alpha, delta, p_norm)
    check_adaptive_parameters_2(; max_attempts)
    total_attempts = MVector{1,Int64}(0)

    return Embedded(epsilon, alpha, delta, p_norm, max_attempts,
                    total_attempts, benchmark_diffeq)
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

get_adaptive_local_order(::CentralDiff, args...) = 2.0

function get_adaptive_local_order(::Doubling,
                                  order::SVector{P,T}) where {P, T <: AbstractFloat}
    return order[1] + 1.0
end

function get_adaptive_local_order(::Embedded,
                                  order::SVector{P,T}) where {P, T <: AbstractFloat}
    return minimum(order) + 1.0
end

function rescale_tolerance(::CentralDiff, order::SVector{P,T}) where {P,T <: AbstractFloat}
    return 1.0 / order[1]
end

function rescale_tolerance(::Doubling, order::SVector{P,T}) where {P, T <: AbstractFloat}
    return order[1] / (1.0 + order[1])
end

function rescale_tolerance(::Embedded, order::SVector{P,T}) where {P, T <: AbstractFloat}
    return minimum(order) / maximum(order)
end

_reconstruct_adaptive(::CentralDiff; kwargs...) = CentralDiff(; kwargs...)

function _reconstruct_adaptive(adaptive::Doubling; kwargs...)
    return Doubling(; kwargs..., max_attempts = adaptive.max_attempts)
end

function _reconstruct_adaptive(adaptive::Embedded; kwargs...)
    return Embedded(; kwargs..., max_attempts = adaptive.max_attempts)
end

function reconstruct_adaptive(adaptive::AdaptiveStepSize, method::ODEMethod,
                              precision::Type{T}) where T <: AbstractFloat

    epsilon = adaptive.epsilon |> precision
    alpha   = adaptive.alpha |> precision
    delta   = adaptive.delta |> precision
    p_norm  = adaptive.p_norm |> precision

    # rescale tolerance parameters
    if !adaptive.benchmark_diffeq
        repower_high = rescale_tolerance(adaptive, method.order)
        epsilon ^= repower_high
        alpha   ^= repower_high
        delta   ^= repower_high
    end

    return _reconstruct_adaptive(adaptive; epsilon, alpha, delta, p_norm)
end

function compute_step_rejection_rate(adaptive::AdaptiveStepSize, timer::TimeLimit)
    if hasproperty(adaptive, :total_attempts)
        @unpack total_steps = timer
        @unpack total_attempts = adaptive
        return 100.0 * (1.0 - total_steps[1]/total_attempts[1])
    end
    return 0.0
end