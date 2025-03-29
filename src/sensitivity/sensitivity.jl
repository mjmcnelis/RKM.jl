
abstract type SensitivityMethod end

struct NoSensitivity <: SensitivityMethod end

@kwdef struct DecoupledDirect{JM, JVM} <: SensitivityMethod where {JM <: JacobianMethod,
                                                                   JVM <: JacobianVectorMethod}
    jacobian_method::JM = FiniteJacobian()
    jacobian_vector_method::JVM = FiniteJacobianVector()
end

function set_jacobian_cache(sensitivity::NoSensitivity, args...)
    return sensitivity
end

# just borrow stage finder for now
function set_jacobian_cache(sensitivity::DecoupledDirect, ode_wrap_p!, f, y, p)
    @unpack jacobian_method = sensitivity
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
    @set! sensitivity.jacobian_method.cache = cache
    return sensitivity
end

function set_jacobian_vector_cache(sensitivity::NoSensitivity, args...)
    return sensitivity
end

function set_jacobian_vector_cache(sensitivity::DecoupledDirect, f, y)
    @unpack jacobian_vector_method = sensitivity

    if jacobian_vector_method isa FiniteJacobianVector
        cache_1 = similar(y)
        cache_2 = similar(y)
        cache_3 = similar(y)
        jvm = FiniteJacobianVector(cache_1, cache_2, cache_3)
    elseif jacobian_vector_method isa ForwardJacobianVector
        cache_1 = Dual{DeivVecTag}.(y, y)
        cache_2 = Dual{DeivVecTag}.(f, f)
        jvm = ForwardJacobianVector(cache_1, cache_2)
    end
    @set! sensitivity.jacobian_vector_method = jvm

    return sensitivity
end

function explicit_sensitivity_stage!(::NoSensitivity, args...)
    return nothing
end

function implicit_sensitivity_stage!(::NoSensitivity, args...)
    return nothing
end

function explicit_sensitivity_stage!(sensitivity, stage_idx, stage_finder, t, dt,
                                     update_cache, ode_wrap!, ode_wrap_p!)

    @unpack y_tmp, f_tmp, J, S_tmp, dS = update_cache
    @unpack jacobian_method, jacobian_vector_method = sensitivity
    # TODO: wrap jvm into a function
    @unpack cache_1, cache_2 = jacobian_vector_method
    if jacobian_vector_method isa FiniteJacobianVector
        # TODO: remove 3rd cache and use f_tmp instead
        @unpack cache_3 = jacobian_vector_method
    end

    # evaluate explicit term dS = dt*(J*S_tmp + df/dp)
    dS_stage = view(dS, :, :, stage_idx)

    ode_wrap!.t[1] = t
    @unpack p = ode_wrap!

    # TODO: wrap this in a function
    # compute state-Jacobian/sensitivity matrix product
    for j in 1:length(p)                # dS <- J*S_tmp
        Jv = view(dS_stage, :, j)
        v = view(S_tmp, :, j)
        # auto_jacvec!(Jv, ode_wrap!, y_tmp, v, cache_1, cache_2)
        num_jacvec_tmp!(Jv, ode_wrap!, y_tmp, v, cache_1, cache_2, cache_3)
        # note: SparseDiffTools version is slow b/c of norm
        # num_jacvec!(Jv, ode_wrap!, y_tmp, v, cache_1, cache_2, cache_3)
    end

    # compute parameter-Jacobian df/dp
    @.. ode_wrap_p!.y = y_tmp
    ode_wrap_p!.t[1] = t
    evaluate_parameter_jacobian!(jacobian_method, S_tmp, ode_wrap_p!, p, f_tmp)

    @.. dS_stage = dt*(dS_stage + S_tmp)# dS <- dt*(J*S_tmp + df/dp)

    return nothing
end

function implicit_sensitivity_stage!(sensitivity, stage_idx, stage_finder, t, dt,
                                     update_cache, ode_wrap!, ode_wrap_p!, A)

    explicit_sensitivity_stage!(sensitivity, stage_idx, stage_finder, t, dt,
                                update_cache, ode_wrap!, ode_wrap_p!)

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

# note: modified version of num_jacvec! from SparseDiffTools
function num_jacvec_tmp!(dy, f, x, v, cache_1, cache_2, cache_3;
                         relstep = sqrt(eps(1.0)), absstep = relstep)
    v_norm = sqrt(sum(abs2, v))

    if v_norm == 0.0
        @.. dy = 0.0
    else
        ϵ = max(absstep, relstep*abs(dot(x, v))/v_norm)
        @.. cache_3 = x + ϵ*v
        f(cache_2, cache_3)
        # note: only needs to be done once (reuse f_tmp?)
        f(cache_1, x)
        @.. dy = (cache_2 - cache_1) / ϵ
    end
    return nothing
end
