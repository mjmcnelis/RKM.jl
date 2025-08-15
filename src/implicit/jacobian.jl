
abstract type JacobianMethod end

@kwdef struct FiniteJacobian{JC, T} <: JacobianMethod where {JC <: JacobianCache,
                                                             T <: AbstractFloat}
    cache::JC = JacobianCache([0.0])
    sparsity::SparseMatrixCSC{T,Int64} = SparseMatrixCSC(Float64[;;])
    evaluations::MVector{1,Int64} = MVector{1,Int64}(0)
end

@kwdef struct ForwardJacobian{JC, T} <: JacobianMethod where {JC <: ForwardColorJacCache,
                                                              T <: AbstractFloat}
    cache::JC = ForwardColorJacCache(nothing, [0.0])
    sparsity::SparseMatrixCSC{T,Int64} = SparseMatrixCSC(Float64[;;])
    evaluations::MVector{1,Int64} = MVector{1,Int64}(0)
end

function reconstruct_jacobian(jacobian_method::FiniteJacobian, ode_wrap!::W,
                              f::Vector{T}, x::Vector{T2}) where {W <: Wrapper,
                                                                  T <: AbstractFloat,
                                                                  T2 <: AbstractFloat}

    evaluations = jacobian_method.evaluations
    sparsity = jacobian_method.sparsity

    evaluations[1] = 0

    if size(sparsity) == (length(f), length(x))
        colorvec = matrix_colors(sparsity)
        cache = JacobianCache(x, f; colorvec, sparsity)
    else
        # TODO: T = Double64, T2 = Float64 doesn't work here
        cache = JacobianCache(x, f)
    end
    @set! jacobian_method.cache = cache

    return jacobian_method
end

function reconstruct_jacobian(jacobian_method::ForwardJacobian, ode_wrap!::W,
                              f::Vector{T}, x::Vector{T2}) where {W <: Wrapper,
                                                                  T <: AbstractFloat,
                                                                  T2 <: AbstractFloat}

    evaluations = jacobian_method.evaluations
    sparsity = jacobian_method.sparsity

    evaluations[1] = 0
    if size(sparsity) == (length(f), length(x))
        colorvec = matrix_colors(sparsity)
        cache = ForwardColorJacCache(ode_wrap!, x; colorvec, sparsity)
    else
        cache = ForwardColorJacCache(ode_wrap!, x)
    end
    @set! jacobian_method.cache = cache

    return jacobian_method
end

function evaluate_jacobian!(jacobian_method::FiniteJacobian,
                            J, ode_wrap!, x, args...)

    cache = jacobian_method.cache
    evaluations = jacobian_method.evaluations

    finite_difference_jacobian!(J, ode_wrap!, x, cache)
    evaluations[1] += 1

    return nothing
end

function evaluate_jacobian!(jacobian_method::ForwardJacobian,
                            J, ode_wrap!, x, args...)

    cache = jacobian_method.cache
    evaluations = jacobian_method.evaluations

    forwarddiff_color_jacobian!(J, ode_wrap!, x, cache)
    evaluations[1] += 1

    return nothing
end

function root_jacobian!(J::Matrix{T}, A::T, dt::T) where T <: AbstractFloat

    @.. J = J * (-A*dt)
    for k in diagind(J)
        J[k] = J[k] + 1.0
    end

    return nothing
end

function root_jacobian!(J::SparseMatrixCSC{T,Int64}, A::T, dt::T) where T <: AbstractFloat

    @.. J.nzval = J.nzval * (-A*dt)
    for k in diagind(J)
        J[k] = J[k] + 1.0
    end

    return nothing
end
