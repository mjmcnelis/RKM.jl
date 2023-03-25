
abstract type Controller end 

struct PIDController{T} <: Controller where {T <: AbstractFloat}
    kI::T
    kP::T
    kD::T
    e_prev::MVector{1,T}
end 

function PIDController(; kI = 0.3, kP = 0.4, kD = 0.0, 
                        precision::Type{T} = Float64) where {T <: AbstractFloat}
    e_prev = MVector{1, precision}(1.0)
    PIDController(kI, kP, kD, e_prev)
end

function rescale_time_step(controller::PIDController, tol::T, 
                           e_norm::T, order::T) where T <: AbstractFloat
    @unpack kI, kP, e_prev = controller 
    kI /= (1.0 + order)
    kP /= (1.0 + order)
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
    return PIDController(; kI, kP, kD, precision)
end