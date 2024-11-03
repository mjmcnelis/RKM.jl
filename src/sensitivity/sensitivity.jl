
abstract type SensitivityMethod end

struct NoSensitivity <: SensitivityMethod end

@kwdef struct DecoupledDirect{JM} <: SensitivityMethod
    jacobian_method::JM = FiniteJacobian()
end

set_jacobian_cache(sensitivity_method::NoSensitivity, args...) = sensitivity_method

# just borrow stage finder for now
function set_jacobian_cache(sensitivity_method::DecoupledDirect, p, y)
    cache = JacobianCache(p, y)
    @set! sensitivity_method.jacobian_method.cache = cache
    return sensitivity_method
end

explicit_sensitivity_stage!(::NoSensitivity, args...) = nothing

function explicit_sensitivity_stage!(sensitivity_method, stage_idx, stage_finder,
                                     t, dt, update_cache, ode_wrap!, ode_wrap_p!, FE
                                    )

    @unpack y_tmp, f_tmp, J, S_tmp, dS = update_cache

    # TOD: put here?
    ode_wrap!.t[1] = t
    ode_wrap_p!.t[1] = t

    @unpack jacobian_method = stage_finder
    cache_p = sensitivity_method.jacobian_method.cache
    @unpack p = ode_wrap!

    # evaluate explicit term dS = dt*(J*S_tmp + df/dp)
    dS_stage = view(dS, :, :, stage_idx)

    # compute Jacobian df/dy
    evaluate_system_jacobian!(jacobian_method, FE, J,
                              ode_wrap!, y_tmp, f_tmp)

    # matrix multiplication is the bottleneck
    mul!(dS_stage, J, S_tmp)            # dS <- J*S_tmp

    # compute Jacobian df/dp
    @.. ode_wrap_p!.y = y_tmp
    finite_difference_jacobian!(S_tmp, ode_wrap_p!, p, cache_p)

    @.. dS_stage = dS_stage + S_tmp     # dS <- dS + df/dp
    @.. dS_stage = dS_stage * dt        # dS <- dS*dt

    return nothing
end