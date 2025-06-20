
abstract type SensitivityMethod end

struct NoSensitivity <: SensitivityMethod end

@kwdef struct DecoupledDirect{JM, JVM} <: SensitivityMethod where {JM <: JacobianMethod,
                                                                   JVM <: JacobianVectorMethod}
    param_jacobian::JM = FiniteJacobian()
    jacobian_vector::JVM = FiniteJacobianVector()
end

function reconstruct_sensitivity(sensitivity::NoSensitivity, args...)
    return sensitivity
end

function reconstruct_sensitivity(sensitivity::DecoupledDirect,
                                 ode_wrap_p!::ODEWrapperParam,
                                 f::Vector{T}, p::Vector{Float64}) where T <: AbstractFloat
    @unpack param_jacobian, jacobian_vector = sensitivity

    param_jacobian = reconstruct_jacobian(param_jacobian, ode_wrap_p!, f, p)
    jacobian_vector = reconstruct_jacobian_vector(jacobian_vector, f)

    @set! sensitivity.param_jacobian = param_jacobian
    @set! sensitivity.jacobian_vector = jacobian_vector
    return sensitivity
end

function explicit_sensitivity_stage!(::NoSensitivity, args...)
    return nothing
end

function explicit_sensitivity_stage!(sensitivity, stage_idx, state_jacobian, t_tmp,
                                     dt, update_cache, ode_wrap_y!, ode_wrap_p!,
                                     method)

    @unpack y_tmp, f_tmp, J, S, S_tmp, dS = update_cache
    @unpack param_jacobian, jacobian_vector = sensitivity

    @unpack A_T = method

    @.. S_tmp = S
    for j in 1:stage_idx-1
        A_T[j,stage_idx] == 0.0 ? continue : nothing
        dS_stage = view(dS,:,:,j)
        @.. S_tmp = S_tmp + A_T[j,stage_idx]*dS_stage
    end

    # set wrapper variables
    set_wrapper!(ode_wrap_y!, t_tmp)
    set_wrapper!(ode_wrap_p!, t_tmp, y_tmp)

    # evaluate explicit term dS = dt*(J*S_tmp + df/dp)
    dS_stage = view(dS, :, :, stage_idx)

    # compute Jacobian-sensitivity product: dS <- J*S_tmp
    evaluate_jacobian_sensitivity!(jacobian_vector, dS_stage, ode_wrap_y!,
                                   state_jacobian, J, S_tmp, y_tmp, f_tmp)

    # compute parameter-Jacobian: S_tmp <- df/dp
    @unpack p = ode_wrap_y!
    # TODO: would it help to store result in a sparse Jp?
    evaluate_jacobian!(param_jacobian, S_tmp, ode_wrap_p!, p, f_tmp)

    # dS <- dt*(J*S_tmp + df/dp)
    @.. dS_stage = dt*(dS_stage + S_tmp)

    return nothing
end

function implicit_sensitivity_stage!(::NoSensitivity, args...)
    return nothing
end

function implicit_sensitivity_stage!(sensitivity, stage_idx, state_jacobian, t_tmp,
                                     dt, update_cache, ode_wrap_y!, ode_wrap_p!, A,
                                     method)

    explicit_sensitivity_stage!(sensitivity, stage_idx, state_jacobian, t_tmp,
                                dt, update_cache, ode_wrap_y!, ode_wrap_p!,
                                method)

    @unpack y_tmp, f_tmp, J, dS = update_cache

    # compute Jacobian df/dy if not already done in explicit sensitivity
    @unpack jacobian_vector = sensitivity
    if !(jacobian_vector isa NaiveJacobianVector)
        evaluate_jacobian!(state_jacobian, J, ode_wrap_y!, y_tmp, f_tmp)
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
