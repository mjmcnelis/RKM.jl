
abstract type SensitivityMethod end

struct NoSensitivity <: SensitivityMethod end

@kwdef struct DecoupledDirect{JM} <: SensitivityMethod
    jacobian_method::JM = FiniteJacobian()
end

set_jacobian_cache(sensitivity_method::NoSensitivity, args...) = sensitivity_method

# just borrow stage finder for now
function set_jacobian_cache(sensitivity_method::DecoupledDirect, ode_wrap_p!, f, y, p)
    @unpack jacobian_method = sensitivity_method
    if jacobian_method isa FiniteJacobian
        cache = JacobianCache(p, y)
    elseif jacobian_method isa ForwardJacobian
        cache = JacobianConfig(ode_wrap_p!, f, p)
    end
    @set! sensitivity_method.jacobian_method.cache = cache
    return sensitivity_method
end

explicit_sensitivity_stage!(::NoSensitivity, args...) = nothing
implicit_sensitivity_stage!(::NoSensitivity, args...) = nothing

function explicit_sensitivity_stage!(sensitivity_method, stage_idx, stage_finder, t,
                                     dt, update_cache, ode_wrap!, ode_wrap_p!, FE)

    @unpack y_tmp, f_tmp, J, S_tmp, dS = update_cache

    # TOD: put here?
    ode_wrap!.t[1] = t
    ode_wrap_p!.t[1] = t

    jacobian_method_y = stage_finder.jacobian_method
    jacobian_method_p = sensitivity_method.jacobian_method
    @unpack p = ode_wrap!

    # evaluate explicit term dS = dt*(J*S_tmp + df/dp)
    dS_stage = view(dS, :, :, stage_idx)

    # compute Jacobian df/dy
    evaluate_system_jacobian!(jacobian_method_y, FE, J,
                              ode_wrap!, y_tmp, f_tmp)

    # matrix multiplication is the bottleneck
    mul!(dS_stage, J, S_tmp)            # dS <- J*S_tmp

    # compute Jacobian df/dp
    @.. ode_wrap_p!.y = y_tmp
    evaluate_parameter_jacobian!(jacobian_method_p, FE, S_tmp,
                                 ode_wrap_p!, p, f_tmp)
    @.. dS_stage = dS_stage + S_tmp     # dS <- dS + df/dp
    @.. dS_stage = dS_stage * dt        # dS <- dS*dt

    return nothing
end

function implicit_sensitivity_stage!(sensitivity_method, stage_idx, stage_finder, t,
                                     dt, update_cache, ode_wrap!, ode_wrap_p!, FE, A)

    explicit_sensitivity_stage!(sensitivity_method, stage_idx, stage_finder, t,
                                dt, update_cache, ode_wrap!, ode_wrap_p!, FE)

    @unpack J, dS = update_cache

    # linear solve for implicit sensitivity stage
    @.. J = J * (-A*dt)                 # J <- I - A.dt.J
    for k in diagind(J)
        J[k] = J[k] + 1.0
    end
    dS_stage = view(dS, :, :, stage_idx)
    F = lu!(J)
    ldiv!(F, dS_stage)                  # dS_stage <- J \ dS_stage

    return nothing
end

function evaluate_parameter_jacobian!(jacobian_method::ForwardJacobian,
                                      FE, S, ode_wrap_p!, p, f)
    @unpack cache#=, evaluations=# = jacobian_method
    jacobian!(S, ode_wrap_p!, f, p, cache)
    FE[1] += ceil(Int64, length(p)/DEFAULT_CHUNK_THRESHOLD)
    # evaluations[1] += 1
    return nothing
end

function evaluate_parameter_jacobian!(jacobian_method::FiniteJacobian,
                                      FE, S, ode_wrap_p!, p, args...)
    @unpack cache#=, evaluations=# = jacobian_method
    finite_difference_jacobian!(S, ode_wrap_p!, p, cache)
    FE[1] += length(p) + 1
    # evaluations[1] += 1
    return nothing
end
