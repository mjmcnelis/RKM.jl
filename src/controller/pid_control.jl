#  TODO: make a TimeStepController struct to wrap PID and other control logic parameters
abstract type Controller end

struct PIDControl{T} <: Controller where {T <: AbstractFloat}
    beta1::T
    beta2::T
    beta3::T
    e_prev::MVector{2,T}
    tol_prev::MVector{2,T}
    dt_prev::MVector{2,T}
    predictive::Bool
    initialized::MVector{1,Bool}
end

function PIDControlBeta(; beta1 = 0.7, beta2 = -0.4, beta3 = 0.0,
                          predictive::Bool = false,
                          precision::Type{T} = Float64) where {T <: AbstractFloat}

    e_prev = MVector{2, precision}(1.0, 1.0)
    tol_prev = MVector{2, precision}(1.0, 1.0)
    dt_prev = MVector{2, precision}(1.0, 1.0)
    initialized = MVector{1, Bool}(false)
    return PIDControl(beta1, beta2, beta3, e_prev, tol_prev, dt_prev, predictive, initialized)
end

function PIDControlK(; kI = 0.3, kP = 0.4, kD = 0.0,
                       predictive::Bool = false,
                       precision::Type{T} = Float64) where {T <: AbstractFloat}
    # relation between beta and k control parameters
    beta1 = kI + kP + kD
    beta2 = -kP - 2kD
    beta3 = kD

    e_prev = MVector{2, precision}(1.0, 1.0)
    tol_prev = MVector{2, precision}(1.0, 1.0)
    dt_prev = MVector{2, precision}(1.0, 1.0)
    initialized = MVector{1, Bool}(false)
    return PIDControl(beta1, beta2, beta3, e_prev, tol_prev, dt_prev, predictive, initialized)
end

function rescale_time_step(controller::PIDControl, tol::T,
                           e_norm::T) where T <: AbstractFloat

    @unpack beta1, beta2, beta3, e_prev, tol_prev, dt_prev, predictive = controller

    # not sure why PID H312 is worse

    # TODO: 2nd factor allocates with Double64
    rescale = (tol/e_norm)^beta1 * (tol_prev[1]/e_prev[1])^beta2 *
              (tol_prev[2]/e_prev[2])^beta3

    if predictive
        # TODO: double check if it's dt/dt_prev[1] or not
        # either way, I probably need a different set of PID parameters
        rescale *= dt_prev[1]/dt_prev[2]
    end
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
    controller.dt_prev[2]  = controller.dt_prev[1]
    controller.dt_prev[1]  = dt
    return nothing
end

function reconstruct_controller(controller::PIDControl,
                                method::ODEMethod, adaptive::AdaptiveStepSize,
                                precision::Type{T}) where T <: AbstractFloat

    @unpack beta1, beta2, beta3, predictive = controller

    local_order = adaptive isa Embedded ? minimum(method.order) + 1.0 :
                  adaptive isa CentralDiff ? 2.0 : method.order[1] + 1.0

    beta1 = precision(beta1 / local_order)
    beta2 = precision(beta2 / local_order)
    beta3 = precision(beta3 / local_order)

    return PIDControlBeta(; beta1, beta2, beta3, predictive, precision)
end