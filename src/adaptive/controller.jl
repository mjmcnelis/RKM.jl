
abstract type Controller end

"""
$(TYPEDEF)

Time step controller based on PID control

# Fields
$(TYPEDFIELDS)
"""
struct TimeStepController{P, T} <: Controller where {P <: PIDControlBeta,
                                                     T <: AbstractFloat}
    pid::P
    e_prev::Vector{T}
    tol_prev::Vector{T}
    dt_prev::Vector{T}
    initialized::MVector{1,Bool}
end

function TimeStepController(precision::Type{T} = Float64;
                            pid::P = PIControl()) where {P <: PIDControlBeta,
                                                         T <: AbstractFloat}
    e_prev   = ones(precision, 2)   # move these to update_cache
    tol_prev = ones(precision, 2)
    dt_prev  = ones(precision, 3)
    initialized = MVector{1, Bool}(false)

    return TimeStepController(pid, e_prev, tol_prev, dt_prev, initialized)
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
                                method::ODEMethod, adaptive::AdaptiveTimeStep,
                                precision::Type{T}) where T <: AbstractFloat
    @unpack pid = controller
    @unpack order = method

    local_order = get_local_order(adaptive, order)
    repower_high = rescale_tolerance(adaptive, order)

    @set! pid.beta1 = pid.beta1 / local_order
    @set! pid.beta2 = pid.beta2 / local_order
    @set! pid.beta3 = pid.beta3 / local_order

    return TimeStepController(precision; pid)
end

function adjust_final_time_steps!(t::Vector{T}, dt::Vector{T},
                                  tf::T) where T <: AbstractFloat
    if dt[2] > tf - t[1]
        dt[2] = tf - t[1]
        dt[1] = dt[2]
    end
    return nothing
end