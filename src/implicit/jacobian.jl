
abstract type JacobianMethod end

struct FiniteJacobian{JC} <: JacobianMethod where JC <: JacobianCache
    sparsity::SparseMatrixCSC{Float64,Int64}
    cache::JC
    evaluations::MVector{1,Int64}
    time_subroutine::Bool
    runtime::MVector{1,Float64}
end

"""
    FiniteJacobian(; sparsity = SparseMatrixCSC(Float64[;;]),)

Outer constructor for `FiniteJacobian`, where you can set a `sparsity` pattern
for the Jacobian matrix.
"""
function FiniteJacobian(; sparsity = SparseMatrixCSC(Float64[;;]),)

    cache = JacobianCache([0.0])
    evaluations = MVector{1,Int64}(0)
    time_subroutine = false
    runtime = MVector{1,Float64}(0.0)

    return FiniteJacobian(sparsity, cache, evaluations, time_subroutine, runtime)
end

struct ForwardJacobian{JC} <: JacobianMethod where JC <: JacobianConfig
    cache::JC
    evaluations::MVector{1,Int64}
    time_subroutine::Bool
    runtime::MVector{1,Float64}
end

"""
    ForwardJacobian()

Outer constructor for `ForwardJacobian`. You cannot pass a `sparsity` pattern
for the Jacobian matrix.
"""
function ForwardJacobian()

    cache = JacobianConfig(nothing, [0.0], [0.0])
    evaluations = MVector{1,Int64}(0)
    time_subroutine = false
    runtime = MVector{1,Float64}(0.0)

    return ForwardJacobian(cache, evaluations, time_subroutine, runtime)
end

struct ForwardColorJacobian{JC} <: JacobianMethod where JC <: ForwardColorJacCache
    sparsity::SparseMatrixCSC{Float64,Int64}
    cache::JC
    evaluations::MVector{1,Int64}
    time_subroutine::Bool
    runtime::MVector{1,Float64}
end

"""
    ForwardColorJacobian(; sparsity = SparseMatrixCSC(Float64[;;]),)

Outer constructor for `ForwardColorJacobian`, where you can set a `sparsity` pattern
for the Jacobian matrix.
"""
function ForwardColorJacobian(; sparsity = SparseMatrixCSC(Float64[;;]),)

    cache = ForwardColorJacCache(nothing, [0.0])
    evaluations = MVector{1,Int64}(0)
    time_subroutine = false
    runtime = MVector{1,Float64}(0.0)

    return ForwardColorJacobian(sparsity, cache, evaluations, time_subroutine, runtime)
end

function Base.show(io::IO, jacobian_method::JM) where JM <: JacobianMethod

    println("$(JM.name.name)")
    println("---------------------")
    if hasproperty(jacobian_method, :sparsity)
        print("sparsity = ")
        display(jacobian_method.sparsity)
    end
    println("cache = $(typeof(jacobian_method.cache))")

    if jacobian_method.evaluations == 0
        println("evaluations = N/A")
    else
        println("evaluations = $(jacobian_method.evaluations)")
    end

    if iszero(jacobian_method.runtime)
        println("runtime = N/A")
    else
        println("runtime = $(jacobian_method.runtime)")
    end

    return nothing
end

function reconstruct_jacobian(jacobian_method::FiniteJacobian, ode_wrap!::W,
                              f::Vector{T}, x::Vector{T},
                              time_subroutine::Bool) where {W <: Wrapper,
                                                            T <: AbstractFloat}
    sparsity = jacobian_method.sparsity
    evaluations = jacobian_method.evaluations
    runtime = jacobian_method.runtime

    evaluations[1] = 0
    runtime[1] = 0.0

    if size(sparsity) == (length(f), length(x))
        colorvec = matrix_colors(sparsity)
        cache = JacobianCache(x, f; colorvec, sparsity)
    else
        cache = JacobianCache(x, f)
    end
    @set! jacobian_method.cache = cache
    @set! jacobian_method.time_subroutine = time_subroutine

    return jacobian_method
end

function reconstruct_jacobian(jacobian_method::ForwardJacobian, ode_wrap!::W,
                              f::Vector{T}, x::Vector{T},
                              time_subroutine::Bool) where {W <: Wrapper,
                                                            T <: AbstractFloat}
    evaluations = jacobian_method.evaluations
    runtime = jacobian_method.runtime

    evaluations[1] = 0
    runtime[1] = 0.0

    cache = JacobianConfig(ode_wrap!, f, x)

    @set! jacobian_method.cache = cache
    @set! jacobian_method.time_subroutine = time_subroutine

    return jacobian_method
end

function reconstruct_jacobian(jacobian_method::ForwardColorJacobian, ode_wrap!::W,
                              f::Vector{T}, x::Vector{T},
                              time_subroutine::Bool) where {W <: Wrapper,
                                                            T <: AbstractFloat}
    sparsity = jacobian_method.sparsity
    evaluations = jacobian_method.evaluations
    runtime = jacobian_method.runtime

    evaluations[1] = 0
    runtime[1] = 0.0

    if size(sparsity) == (length(f), length(x))
        colorvec = matrix_colors(sparsity)
        cache = ForwardColorJacCache(ode_wrap!, x; colorvec, sparsity)
    else
        # note: doesn't work if f and x have different dimensions
        cache = ForwardColorJacCache(ode_wrap!, x)
    end
    @set! jacobian_method.cache = cache
    @set! jacobian_method.time_subroutine = time_subroutine

    return jacobian_method
end

function evaluate_jacobian!(jacobian_method::FiniteJacobian,
                            J, ode_wrap!, x, args...)

    cache = jacobian_method.cache
    evaluations = jacobian_method.evaluations
    time_subroutine = jacobian_method.time_subroutine
    runtime = jacobian_method.runtime

    if time_subroutine && evaluations[1] % SAMPLE_INTERVAL == 0
        stats = @timed finite_difference_jacobian!(J, ode_wrap!, x, cache)
        runtime[1] += SAMPLE_INTERVAL*stats.time
    else
        finite_difference_jacobian!(J, ode_wrap!, x, cache)
    end
    evaluations[1] += 1

    return nothing
end

function evaluate_jacobian!(jacobian_method::ForwardJacobian,
                            J, ode_wrap!, x, f)

    cache = jacobian_method.cache
    evaluations = jacobian_method.evaluations

    time_subroutine = jacobian_method.time_subroutine
    runtime = jacobian_method.runtime

    if time_subroutine && evaluations[1] % SAMPLE_INTERVAL == 0
        stats = @timed jacobian!(J, ode_wrap!, f, x, cache)
        runtime[1] += SAMPLE_INTERVAL*stats.time
    else
        jacobian!(J, ode_wrap!, f, x, cache)
    end
    evaluations[1] += 1

    return nothing
end

function evaluate_jacobian!(jacobian_method::ForwardColorJacobian,
                            J, ode_wrap!, x, args...)

    cache = jacobian_method.cache
    evaluations = jacobian_method.evaluations

    time_subroutine = jacobian_method.time_subroutine
    runtime = jacobian_method.runtime

    if time_subroutine && evaluations[1] % SAMPLE_INTERVAL == 0
        stats = @timed forwarddiff_color_jacobian!(J, ode_wrap!, x, cache)
        runtime[1] += SAMPLE_INTERVAL*stats.time
    else
        forwarddiff_color_jacobian!(J, ode_wrap!, x, cache)
    end
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
