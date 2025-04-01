
abstract type SensitivityMethod end

struct NoSensitivity <: SensitivityMethod end

@kwdef struct DecoupledDirect{JM, JVM} <: SensitivityMethod where {JM <: JacobianMethod,
                                                                   JVM <: JacobianVectorMethod}
    param_jacobian::JM = FiniteJacobian()
    jacobian_vector::JVM = FiniteJacobianVector()
end

function set_jacobian_cache(sensitivity::NoSensitivity, args...)
    return sensitivity
end

# just borrow stage finder for now
function set_jacobian_cache(sensitivity::DecoupledDirect, ode_wrap_p!::ODEWrapperParam,
                            f::Vector{T}, y::Vector{T},
                            p::Vector{Float64}) where T <: AbstractFloat

    @unpack param_jacobian = sensitivity
    if param_jacobian isa FiniteJacobian
        @unpack sparsity = param_jacobian
        if size(sparsity) == (length(y), length(p))
            colorvec = matrix_colors(sparsity)
            cache = JacobianCache(p, y; colorvec, sparsity)
        else
            cache = JacobianCache(p, y)
        end
    elseif param_jacobian isa ForwardJacobian
        cache = JacobianConfig(ode_wrap_p!, f, p)
    elseif param_jacobian isa ForwardColorJacobian
        @unpack sparsity = param_jacobian
        if size(sparsity) == (length(y), length(p))
            colorvec = matrix_colors(sparsity)
            cache = ForwardColorJacCache(ode_wrap_p!, p; colorvec, sparsity)
        else
            cache = ForwardColorJacCache(ode_wrap_p!, p)
        end
    end
    @set! sensitivity.param_jacobian.cache = cache
    return sensitivity
end

function set_jacobian_vector_cache(sensitivity::NoSensitivity, args...)
    return sensitivity
end

function reconstruct_jacobian_vector(::FiniteJacobianVector, y, args...)
    cache_1 = similar(y)
    cache_2 = similar(y)
    return FiniteJacobianVector(cache_1, cache_2)
end

function reconstruct_jacobian_vector(::ForwardJacobianVector, y, f)
    cache_1 = Dual{DeivVecTag}.(y, y)
    cache_2 = Dual{DeivVecTag}.(f, f)
    return ForwardJacobianVector(cache_1, cache_2)
end

# TODO: replace w/ main reconstruction where remake both J and Jv
function set_jacobian_vector_cache(sensitivity::DecoupledDirect, y, f)
    @unpack jacobian_vector = sensitivity

    if jacobian_vector isa NaiveJacobianVector
        return sensitivity
    end
    jvm = reconstruct_jacobian_vector(jacobian_vector, y, f)
    @set! sensitivity.jacobian_vector = jvm

    return sensitivity
end

function explicit_sensitivity_stage!(::NoSensitivity, args...)
    return nothing
end

function explicit_sensitivity_stage!(sensitivity, stage_idx, stage_finder, t, dt,
                                     update_cache, ode_wrap!, ode_wrap_p!)

    @unpack y_tmp, f_tmp, J, S_tmp, dS = update_cache
    @unpack param_jacobian, jacobian_vector = sensitivity

    # evaluate explicit term dS = dt*(J*S_tmp + df/dp)
    dS_stage = view(dS, :, :, stage_idx)

    # compute Jacobian-sensitivity product: dS <- J*S_tmp
    # TODO: maybe just pass update_cache
    evaluate_jacobian_sensitivity!(jacobian_vector, dS_stage, ode_wrap!,
                                   stage_finder, t, J, S_tmp, y_tmp, f_tmp)

    # compute parameter-Jacobian: S_tmp <- df/dp
    @unpack p = ode_wrap!
    # TODO: wrapper utils function to set variables could be useful
    ode_wrap_p!.t[1] = t
    @.. ode_wrap_p!.y = y_tmp
    # TODO: would it help to store result in a sparse Jp?
    evaluate_parameter_jacobian!(param_jacobian, S_tmp, ode_wrap_p!, p, f_tmp)

    # dS <- dt*(J*S_tmp + df/dp)
    @.. dS_stage = dt*(dS_stage + S_tmp)

    return nothing
end

function implicit_sensitivity_stage!(::NoSensitivity, args...)
    return nothing
end

function implicit_sensitivity_stage!(sensitivity, stage_idx, stage_finder, t, dt,
                                     update_cache, ode_wrap!, ode_wrap_p!, A)

    explicit_sensitivity_stage!(sensitivity, stage_idx, stage_finder, t, dt,
                                update_cache, ode_wrap!, ode_wrap_p!)

    @unpack y_tmp, f_tmp, J, dS = update_cache

    # compute Jacobian df/dy if not already done in explicit sensitivity
    @unpack jacobian_vector = sensitivity
    if !(jacobian_vector isa NaiveJacobianVector)
        @unpack jacobian_method = stage_finder
        ode_wrap!.t[1] = t
        evaluate_system_jacobian!(jacobian_method, J, ode_wrap!, y_tmp, f_tmp)
    end

    # linear solve for implicit sensitivity stage
    if J isa SparseMatrixCSC            # J <- I - A.dt.J
        @.. J.nzval = J.nzval * (-A*dt)
    else
        @.. J = J * (-A*dt)
    end
    for k in diagind(J)
        J[k] = J[k] + 1.0
    end
    dS_stage = view(dS, :, :, stage_idx)

    # note: use LinearAlgebra (LinearSolve doesn't support matrix-matrix equations AX = B)
    if J isa SparseMatrixCSC          # dS_stage <- J \ dS_stage
        dS_stage = J \ dS_stage
    else
        F = lu!(J)
        ldiv!(F, dS_stage)
    end

    return nothing
end
