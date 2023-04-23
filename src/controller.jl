#  TODO: make a TimeStepController struct to wrap PID and other control logic parameters
abstract type Controller end

# TODO: just rename as PID
struct PIDController{T} <: Controller where {T <: AbstractFloat}
    kI::T
    kP::T
    kD::T
    e_prev::MVector{2,T}
    dt_prev::MVector{2,T}
    predictive::Bool
end 

function PIDControllerBeta(; beta1 = 0.7, beta2 = -0.4, beta3 = 0.0, 
                             predictive::Bool = false,
                             precision::Type{T} = Float64) where {T <: AbstractFloat}
    # relation between k and beta control parameters
    kI = beta1 + beta2 + beta3
    kP = -beta2 - 2beta3
    kD = beta3

    e_prev = MVector{2, precision}(1.0, 1.0)
    dt_prev = MVector{2, precision}(1.0, 1.0)
    return PIDController(kI, kP, kD, e_prev, dt_prev, predictive)
end

function PIDControllerK(; kI = 0.3, kP = 0.4, kD = 0.0,
                          predictive::Bool = false,
                          precision::Type{T} = Float64) where {T <: AbstractFloat}
    e_prev = MVector{2, precision}(1.0, 1.0)
    dt_prev = MVector{2, precision}(1.0, 1.0)
    return PIDController(kI, kP, kD, e_prev, dt_prev, predictive)
end

function rescale_time_step(controller::PIDController, tol::T, 
                           e_norm::T, order::T) where T <: AbstractFloat
    @unpack kI, kP, kD, e_prev, dt_prev, predictive = controller 
    kI /= (1.0 + order)
    kP /= (1.0 + order)
    kD /= (1.0 + order)

    # @show  dt/dt_prev[1]
    # not sure why PID H312 is worse, maybe I set previous variables wrong?
    # or maybe tolerance should also grab previous values

    # TODO: 2nd factor allocates with Double64 
    rescale = (tol/e_norm)^kI * (e_prev[1]/e_norm)^kP * (e_prev[1]^2/e_prev[2]/e_norm)^kD
    if predictive
        # TODO: double check if it's dt/dt_prev[1] or not 
        # either way, I probably need a different set of PID parameters
        rescale *= dt_prev[1]/dt_prev[2]
    end
    return rescale
end

function initialize_controller!(controller::PIDController, 
                                e_norm::T, dt::T) where T <: AbstractFloat
    controller.e_prev  .= e_norm
    controller.dt_prev .= dt 
    return nothing 
end

function set_previous_control_vars!(controller::PIDController, 
                                    e_norm::T, dt::T) where T <: AbstractFloat
    controller.e_prev[2]  = controller.e_prev[1]
    controller.e_prev[1]  = e_norm
    controller.dt_prev[2] = controller.dt_prev[1]
    controller.dt_prev[1] = dt
    return nothing 
end

function reconstruct_controller(controller::PIDController, 
                                precision::Type{T}) where T <: AbstractFloat
    @unpack kI, kP, kD, predictive = controller
    kI = precision(kI) 
    kP = precision(kP) 
    kD = precision(kD) 
    return PIDControllerK(; kI, kP, kD, predictive, precision)
end