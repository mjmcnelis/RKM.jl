
abstract type Controller end

struct PIDController{T} <: Controller where {T <: AbstractFloat}
    kI::T
    kP::T
    kD::T
    e_prev::MVector{2,T}
end 

function PIDControllerBeta(; beta1 = 0.7, beta2 = -0.4, beta3 = 0.0,
                             precision::Type{T} = Float64) where {T <: AbstractFloat}
    # relation between k and beta control parameters
    kI = beta1 + beta2 + beta3
    kP = -beta2 - 2beta3
    kD = beta3
    e_prev = MVector{2, precision}(1.0, 1.0)
    return PIDController(kI, kP, kD, e_prev)
end

function PIDControllerK(; kI = 0.3, kP = 0.4, kD = 0.0,
                          precision::Type{T} = Float64) where {T <: AbstractFloat}
    e_prev = MVector{2, precision}(1.0, 1.0)
    return PIDController(kI, kP, kD, e_prev)
end

function rescale_time_step(controller::PIDController, tol::T, 
                           e_norm::T, order::T) where T <: AbstractFloat
    @unpack kI, kP, e_prev = controller 
    kI /= (1.0 + order)
    kP /= (1.0 + order)
    # TODO: add kD controller term (need formula worked out) 
    # also need to set previous-previous error e_{n-2}

    # TODO: 2nd factor allocates with Double64 
    rescale = (tol / e_norm)^kI * (e_prev[1] / e_norm)^kP
    return rescale
end

function set_previous_control_error!(controller::PIDController, 
                                     e_norm::T) where T <: AbstractFloat
    controller.e_prev[1] = e_norm
    return nothing 
end

function reconstruct_controller(controller::PIDController, 
                                precision::Type{T}) where T <: AbstractFloat
    @unpack kI, kP, kD = controller
    kI = precision(kI) 
    kP = precision(kP) 
    kD = precision(kD) 
    return PIDControllerK(; kI, kP, kD, precision)
end