
abstract type Controller end

"""
$(TYPEDEF)

Time step controller based on PID control

# Fields
$(TYPEDFIELDS)
"""
struct TimeStepController{P, L, T} <: Controller where {P <: PIDControlBeta,
                                                        L <: LimiterMethod,
                                                        T <: AbstractFloat}
    pid::P
    limiter::L
    e_prev::VectorMVector{2,T}
    tol_prev::VectorMVector{2,T}
    dt_prev::VectorMVector{3,T}
    initialized::MVector{1,Bool}
end

function TimeStepController(; pid = PIControl(), limiter = PiecewiseLimiter(),
                              precision::Type{T} = Float64) where {T <: AbstractFloat}
    # initialize vectors
    if precision == BigFloat
        e_prev   = ones(precision, 2)
        tol_prev = ones(precision, 2)
        dt_prev  = ones(precision, 3)
    else
        e_prev   = MVector{2, precision}(1.0, 1.0)
        tol_prev = MVector{2, precision}(1.0, 1.0)
        dt_prev  = MVector{3, precision}(1.0, 1.0, 1.0)
    end

    initialized = MVector{1, Bool}(false)

    return TimeStepController(pid, limiter, e_prev, tol_prev, dt_prev, initialized)
end

function rescale_time_step(controller::TimeStepController, tol::T,
                           e_norm::T) where T <: AbstractFloat

    @unpack pid, e_prev, tol_prev, dt_prev = controller
    @unpack beta1, beta2, beta3, alpha2, alpha3 = pid

    # note: use a slightly modified ^ operator for DoubleFloats (see explog.jl)
    rescale = (tol/e_norm)^beta1 * (tol_prev[1]/e_prev[1])^beta2 *
              (tol_prev[2]/e_prev[2])^beta3 * (dt_prev[2]/dt_prev[1])^alpha2 *
              (dt_prev[3]/dt_prev[2])^alpha3

    return rescale
end

function initialize_controller!(controller::TimeStepController, e_norm::T,
                                tol::T, dt::T) where T <: AbstractFloat
    controller.e_prev   .= e_norm
    controller.tol_prev .= tol
    controller.dt_prev  .= dt
    return nothing
end

function set_previous_control_vars!(controller::TimeStepController, e_norm::T,
                                    tol::T, dt::T) where T <: AbstractFloat
    controller.e_prev[2]   = controller.e_prev[1]
    controller.e_prev[1]   = e_norm
    controller.tol_prev[2] = controller.tol_prev[1]
    controller.tol_prev[1] = tol
    controller.dt_prev[3]  = controller.dt_prev[2]
    controller.dt_prev[2]  = controller.dt_prev[1]
    controller.dt_prev[1]  = dt
    return nothing
end

function reconstruct_controller(controller::TimeStepController,
                                method::ODEMethod, adaptive::AdaptiveStepSize,
                                precision::Type{T}) where T <: AbstractFloat

    @unpack pid, limiter = controller
    @unpack beta1, beta2, beta3, alpha2, alpha3 = pid
    @unpack safety, high = limiter
    @unpack order = method

    local_order = adaptive isa Embedded ? minimum(order) + 1.0 :
                  adaptive isa CentralDiff ? 2.0 :
                  adaptive isa Doubling ? order[1] + 1.0 : error()

    repower_high = adaptive isa Embedded ? minimum(order)/maximum(order) :
                   adaptive isa CentralDiff ? 1.0/order[1] :
                   adaptive isa Doubling ? order[1]/(1.0 + order[1]) : error()

    @set! pid.beta1 = beta1 / local_order
    @set! pid.beta2 = beta2 / local_order
    @set! pid.beta3 = beta3 / local_order
    @set! limiter.high = high^repower_high

    # TODO: add message
    # reminder: reason I need this is b/c I default rescale = high when error = 0
    @assert safety * limiter.high > 1.0

    return TimeStepController(; pid, limiter, precision)
end