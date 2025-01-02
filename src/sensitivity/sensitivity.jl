
abstract type SensitivityMethod end

struct NoSensitivity <: SensitivityMethod end

@kwdef struct DecoupledDirect{JM} <: SensitivityMethod
    jacobian_method::JM = FiniteJacobian()
end

function set_jacobian_cache(sensitivity_method::NoSensitivity, args...)
    return sensitivity_method
end

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

function explicit_sensitivity_stage!(::NoSensitivity, args...)
    return nothing
end

function implicit_sensitivity_stage!(::NoSensitivity, args...)
    return nothing
end

function explicit_sensitivity_stage!(sensitivity_method, stage_idx, stage_finder, t,
                                     dt, update_cache, ode_wrap!, ode_wrap_p!)

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
    evaluate_system_jacobian!(jacobian_method_y, J, ode_wrap!, y_tmp, f_tmp)

    # matrix multiplication is the bottleneck
    mul!(dS_stage, J, S_tmp)            # dS <- J*S_tmp

    # compute Jacobian df/dp
    @.. ode_wrap_p!.y = y_tmp
    evaluate_parameter_jacobian!(jacobian_method_p, S_tmp, ode_wrap_p!, p, f_tmp)

    @.. dS_stage = dS_stage + S_tmp     # dS <- dS + df/dp
    @.. dS_stage = dS_stage * dt        # dS <- dS*dt

    return nothing
end

function implicit_sensitivity_stage!(sensitivity_method, stage_idx, stage_finder, t,
                                     dt, update_cache, ode_wrap!, ode_wrap_p!, A)

    explicit_sensitivity_stage!(sensitivity_method, stage_idx, stage_finder, t,
                                dt, update_cache, ode_wrap!, ode_wrap_p!)

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
                                      S, ode_wrap_p!, p, f)
    @unpack cache#=, evaluations=# = jacobian_method
    jacobian!(S, ode_wrap_p!, f, p, cache)
    # evaluations[1] += 1
    return nothing
end

function evaluate_parameter_jacobian!(jacobian_method::FiniteJacobian,
                                      S, ode_wrap_p!, p, args...)
    @unpack cache#=, evaluations=# = jacobian_method
    finite_difference_jacobian!(S, ode_wrap_p!, p, cache)
    # evaluations[1] += 1
    return nothing
end
