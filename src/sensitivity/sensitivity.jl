
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

    param_jacobian = reconstruct_jacobian_method(param_jacobian, ode_wrap_p!, f, p)
    jacobian_vector = reconstruct_jacobian_vector(jacobian_vector, f)

    @set! sensitivity.param_jacobian = param_jacobian
    @set! sensitivity.jacobian_vector = jacobian_vector
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
    evaluate_jacobian!(param_jacobian, S_tmp, ode_wrap_p!, p, f_tmp)

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
        @unpack state_jacobian = stage_finder
        ode_wrap!.t[1] = t
        evaluate_jacobian!(state_jacobian, J, ode_wrap!, y_tmp, f_tmp)
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
