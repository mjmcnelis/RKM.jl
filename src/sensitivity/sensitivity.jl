
abstract type SensitivityMethod end

struct NoSensitivity <: SensitivityMethod end

@kwdef struct DecoupledDirect{JM, JVM} <: SensitivityMethod where {JM <: JacobianMethod,
                                                                   JVM <: JacobianVectorMethod}
    jacobian_method::JM = FiniteJacobian()
    # TODO: default to FiniteJacobianVector once have method
    jacobian_vector_method::JVM = ForwardJacobianVector()
end

function set_jacobian_cache(sensitivity_method::NoSensitivity, args...)
    return sensitivity_method
end

# just borrow stage finder for now
function set_jacobian_cache(sensitivity_method::DecoupledDirect, ode_wrap_p!, f, y, p)
    @unpack jacobian_method = sensitivity_method
    if jacobian_method isa FiniteJacobian
        @unpack sparsity = jacobian_method
        if size(sparsity) == (length(y), length(p))
            colorvec = matrix_colors(sparsity)
            cache = JacobianCache(p, y; colorvec, sparsity)
        else
            cache = JacobianCache(p, y)
        end
    elseif jacobian_method isa ForwardJacobian
        cache = JacobianConfig(ode_wrap_p!, f, p)
    end
    @set! sensitivity_method.jacobian_method.cache = cache
    return sensitivity_method
end

function set_jacobian_vector_cache(sensitivity_method::NoSensitivity, args...)
    return sensitivity_method
end

function set_jacobian_vector_cache(sensitivity_method::DecoupledDirect, f, y)
    @unpack jacobian_vector_method = sensitivity_method

    if jacobian_vector_method isa ForwardJacobianVector
        cache_1 = Dual{DeivVecTag}.(y, y)
        cache_2 = Dual{DeivVecTag}.(f, f)
        @set! sensitivity_method.jacobian_vector_method.cache_1 = cache_1
        @set! sensitivity_method.jacobian_vector_method.cache_2 = cache_2
    end

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
    @unpack jacobian_method, jacobian_vector_method = sensitivity_method
    @unpack cache_1, cache_2 = jacobian_vector_method

    # evaluate explicit term dS = dt*(J*S_tmp + df/dp)
    dS_stage = view(dS, :, :, stage_idx)

    # v = view(S_tmp, :, 1)
    # cache_1 = similar(v)
    # cache_2 = similar(v)
    # cache_3 = similar(v)

    ode_wrap!.t[1] = t
    @unpack p = ode_wrap!

    # TODO: wrap this in a function
    # compute state-Jacobian/sensitivity matrix product
    for j in 1:length(p)                # dS <- J*S_tmp
        Jv = view(dS_stage, :, j)
        v = view(S_tmp, :, j)
        auto_jacvec!(Jv, ode_wrap!, y_tmp, v, cache_1, cache_2)

        # TODO: finite-diff version needs work
        # num_jacvec!(Jv, ode_wrap!, y_tmp, v, cache_1, cache_2, cache_3)
        #=
        # doesn't give me exactly what I want
        finite_difference_jacobian!(Jv, (g,λ) -> (ode_wrap!(g, y_tmp .+ λ.*v)),
                                    zeros(length(y_tmp)))
        =#
    end

    # compute parameter-Jacobian df/dp
    @.. ode_wrap_p!.y = y_tmp
    ode_wrap_p!.t[1] = t
    evaluate_parameter_jacobian!(jacobian_method, S_tmp, ode_wrap_p!, p, f_tmp)

    @.. dS_stage = dt*(dS_stage + S_tmp)# dS <- dt*(J*S_tmp + df/dp)

    return nothing
end

function implicit_sensitivity_stage!(sensitivity_method, stage_idx, stage_finder, t,
                                     dt, update_cache, ode_wrap!, ode_wrap_p!, A)

    explicit_sensitivity_stage!(sensitivity_method, stage_idx, stage_finder, t,
                                dt, update_cache, ode_wrap!, ode_wrap_p!)

    @unpack y_tmp, f_tmp, J, dS = update_cache

    # compute Jacobian df/dy
    @unpack jacobian_method = stage_finder
    ode_wrap!.t[1] = t
    evaluate_system_jacobian!(jacobian_method, J, ode_wrap!, y_tmp, f_tmp)

    # linear solve for implicit sensitivity stage
    @.. J = J * (-A*dt)                 # J <- I - A.dt.J
    for k in diagind(J)
        J[k] = J[k] + 1.0
    end
    dS_stage = view(dS, :, :, stage_idx)
    # note: factorization very slow if J is sparse
    F = lu!(J)
    ldiv!(F, dS_stage)                  # dS_stage <- J \ dS_stage
    # dS_stage = J \ dS_stage

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
