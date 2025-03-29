
abstract type JacobianMethod end

@kwdef struct ForwardJacobian{JC} <: JacobianMethod where JC <: JacobianConfig
    cache::JC = JacobianConfig(nothing, [0.0], [0.0])
    evaluations::MVector{1,Int64} = MVector{1,Int64}(0)
end

@kwdef struct ForwardColorJacobian{JC, T} <: JacobianMethod where {JC <: ForwardColorJacCache,
                                                                   T <: AbstractFloat}
    cache::JC = ForwardColorJacCache(nothing, [0.0])
    sparsity::SparseMatrixCSC{T,Int64} = SparseMatrixCSC(Float64[;;])
    evaluations::MVector{1,Int64} = MVector{1,Int64}(0)
end

@kwdef struct FiniteJacobian{JC, T} <: JacobianMethod where {JC <: JacobianCache,
                                                             T <: AbstractFloat}
    cache::JC = JacobianCache([0.0])
    sparsity::SparseMatrixCSC{T,Int64} = SparseMatrixCSC(Float64[;;])
    evaluations::MVector{1,Int64} = MVector{1,Int64}(0)
end

function evaluate_system_jacobian!(jacobian_method::ForwardJacobian,
                                   J, ode_wrap!, y, f)
    @unpack cache, evaluations = jacobian_method
    jacobian!(J, ode_wrap!, f, y, cache)
    evaluations[1] += 1
    return nothing
end

function evaluate_system_jacobian!(jacobian_method::ForwardColorJacobian,
                                   J, ode_wrap!, y, args...)
    @unpack cache, evaluations = jacobian_method
    forwarddiff_color_jacobian!(J, ode_wrap!, y, cache)
    evaluations[1] += 1
    return nothing
end

function evaluate_system_jacobian!(jacobian_method::FiniteJacobian,
                                   J, ode_wrap!, y, args...)
    @unpack cache, evaluations = jacobian_method
    finite_difference_jacobian!(J, ode_wrap!, y, cache)
    evaluations[1] += 1
    return nothing
end

function evaluate_parameter_jacobian!(param_jacobian::ForwardJacobian,
                                      S, ode_wrap_p!, p, f)
    @unpack cache, evaluations = param_jacobian
    jacobian!(S, ode_wrap_p!, f, p, cache)
    evaluations[1] += 1
    return nothing
end

function evaluate_parameter_jacobian!(param_jacobian::FiniteJacobian,
                                      S, ode_wrap_p!, p, args...)
    @unpack cache, evaluations = param_jacobian
    finite_difference_jacobian!(S, ode_wrap_p!, p, cache)
    evaluations[1] += 1
    return nothing
end
