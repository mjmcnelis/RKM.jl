
abstract type Controller end

"""
$(TYPEDEF)

Time step controller based on PID control

# Fields
$(TYPEDFIELDS)
"""
struct PIDControl{T} <: Controller where {T <: AbstractFloat}
    # why not just always use Float64 for betas?
    beta1::T
    beta2::T
    beta3::T
    alpha2::T
    alpha3::T
    # TODO: wrap safety, low, high in limiter method
    """Safety factor to scale down estimate for predicted time step"""
    safety::T
    """Lower bound on the time step's rate of change"""
    low::T
    """Upper bound on the time step's rate of change"""
    high::T                 # note: CentralDiff used high = 1.5 instead of 5
    e_prev::MVector{2,T}
    tol_prev::MVector{2,T}
    dt_prev::MVector{3,T}
    initialized::MVector{1,Bool}
end

function PIDControlBeta(; beta1 = 0.7, beta2 = -0.4, beta3 = 0.0,
                          alpha2 = 0.0, alpha3 = 0.0,
                          # TODO: come up with limiter method struct to hold safety, low, high
                          safety = 0.8, low = 0.2, high = 5.0,
                          precision::Type{T} = Float64) where {T <: AbstractFloat}

    e_prev = MVector{2, precision}(1.0, 1.0)
    tol_prev = MVector{2, precision}(1.0, 1.0)
    dt_prev = MVector{3, precision}(1.0, 1.0, 1.0)
    initialized = MVector{1, Bool}(false)

    # TODO: once have limiter struct, wrap asserts in function
    @assert 0.0 <= low < 1.0 "low = $low is out of bounds [0, 1)"
    @assert 1.0 < high <= Inf "high = $high is out of bounds (1, Inf]"
    @assert 0.0 < safety < 1.0 "safety = $safety is out of bounds (0, 1)"
    @assert safety*high > 1.0 "safety*high = $(safety*high) is not greater than 1"

    return PIDControl(beta1, beta2, beta3, alpha2, alpha3,
                      safety, low, high,
                      e_prev, tol_prev, dt_prev, initialized)
end

function PIDControlK(; kI = 0.3, kP = 0.4, kD = 0.0,
                       alpha2 = 0.0, alpha3 = 0.0,
                       safety = 0.8, low = 0.2, high = 5.0,
                       precision::Type{T} = Float64) where {T <: AbstractFloat}
    # relation between beta and k control parameters
    beta1 = kI + kP + kD
    beta2 = -kP - 2kD
    beta3 = kD

    return PIDControlBeta(; beta1, beta2, beta3, alpha2, alpha3,
                            safety, low, high, precision)
end

function rescale_time_step(controller::PIDControl, tol::T,
                           e_norm::T) where T <: AbstractFloat

    @unpack beta1, beta2, beta3, alpha2, alpha3, e_prev, tol_prev, dt_prev = controller

    # TODO: 2nd factor allocates with Double64
    rescale = (tol/e_norm)^beta1 * (tol_prev[1]/e_prev[1])^beta2 *
              (tol_prev[2]/e_prev[2])^beta3 * (dt_prev[2]/dt_prev[1])^alpha2 *
              (dt_prev[3]/dt_prev[2])^alpha3

    return rescale
end

function initialize_controller!(controller::PIDControl, e_norm::T,
                                tol::T, dt::T) where T <: AbstractFloat
    controller.e_prev   .= e_norm
    controller.tol_prev .= tol
    controller.dt_prev  .= dt
    return nothing
end

function set_previous_control_vars!(controller::PIDControl, e_norm::T,
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

function reconstruct_controller(controller::PIDControl,
                                method::ODEMethod, adaptive::AdaptiveStepSize,
                                precision::Type{T}) where T <: AbstractFloat

    @unpack beta1, beta2, beta3, alpha2, alpha3, safety, low, high = controller
    @unpack order = method

    local_order = adaptive isa Embedded ? minimum(order) + 1.0 :
                  adaptive isa CentralDiff ? 2.0 :
                  adaptive isa Doubling ? order[1] + 1.0 : error()

    repower_high = adaptive isa Embedded ? minimum(order)/maximum(order) :
                   adaptive isa CentralDiff ? 1.0/order[1] :
                   adaptive isa Doubling ? order[1]/(1.0 + order[1]) : error()

    beta1 = precision(beta1 / local_order)
    beta2 = precision(beta2 / local_order)
    beta3 = precision(beta3 / local_order)
    alpha2 = precision(alpha2)
    alpha3 = precision(alpha3)
    safety = precision(safety)
    low    = precision(low)
    high   = precision(high^repower_high)

    return PIDControlBeta(; beta1, beta2, beta3, alpha2, alpha3,
                            safety, low, high, precision)
end