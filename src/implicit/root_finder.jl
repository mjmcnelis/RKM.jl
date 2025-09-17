
abstract type RootFinderMethod end

@kwdef struct FixedPoint <: RootFinderMethod
    epsilon::Float64 = 1e-8
    p_norm::Float64 = 2.0
    max_iterations::Int64 = 10
    # TODO: track how many iterations were performed each step

    function FixedPoint(epsilon, p_norm, max_iterations)
        @assert p_norm >= 1 "p_norm = $p_norm is not valid"
        return new(epsilon, p_norm, max_iterations)
    end
end

struct Newton{LC} <: RootFinderMethod where LC <: LinearCache
    linear_cache::LC
    epsilon::Float64
    p_norm::Float64
    max_iterations::Int64
    # TODO: track how many iterations were performed each step
    #       may want to have option to bins statistics
    evaluations::MVector{1,Int64}   # note: counter only for runtime sampling
    time_subroutine::Bool
    runtime::MVector{1,Float64}
end

function Newton(; linear_method::AF = LUFactorization(),
                  epsilon::Float64 = 1e-8, p_norm::Float64 = 2.0,
                  max_iterations::Int64 = 10) where AF <: AbstractFactorization

    @assert p_norm >= 1 "p_norm = $p_norm is not valid"

    # configure dummy linear cache (see src/common.jl in LinearSolve.jl)
    J = zeros(0, 0)
    res = zeros(0)
    linear_cache = init(LinearProblem(J, res), linear_method,
                        alias = LinearAliasSpecifier(; alias_A = true, alias_b = true))
    # dummy defaults
    evaluations = MVector{1,Int64}(0)
    time_subroutine = false
    runtime = MVector{1,Float64}(0.0)

    return Newton(linear_cache, epsilon, p_norm, max_iterations,
                  evaluations, time_subroutine, runtime)
end

function reconstruct_root_finder(root_finder::FixedPoint, args...)
    return root_finder
end

function reconstruct_root_finder(root_finder::Newton, res::Vector{T},
                                 J::Union{Matrix{T}, SparseMatrixCSC{T,Int64}},
                                 time_subroutine::Bool) where T <: AbstractFloat

    evaluations = root_finder.evaluations
    runtime = root_finder.runtime

    evaluations[1] = 0
    runtime[1] = 0.0

    linear_method = root_finder.linear_cache.alg
    warn_linear_solver(J, linear_method)

    linear_cache = init(LinearProblem(J, res), linear_method;
                        alias = LinearAliasSpecifier(; alias_A = true, alias_b = true),)

    @set! root_finder.linear_cache = linear_cache
    @set! root_finder.time_subroutine = time_subroutine

    return root_finder
end

function warn_linear_solver(J::Union{Matrix{T}, SparseMatrixCSC{T,Int64}},
                            linear_method::AF) where {T <: AbstractFloat,
                                                      AF <: AbstractFactorization}

    if J isa Matrix && linear_method isa AbstractSparseFactorization
        @warn "Jacobian matrix is dense and may not be compatible with \
               the sparse linear solver $(typeof(linear_method)), either \
               pass a sparsity pattern or use a dense linear solver"
    end
    if T != Float64
        if J isa SparseMatrixCSC && linear_method isa AbstractSparseFactorization
            @warn "Sparse linear solver $(typeof(linear_method)) may not \
                   have support for precision = $T, use Float64 instead"
        else
            @warn "Linear solver $(typeof(linear_method)) may not be \
                   optimized for precision = $T, use Float64 instead"
        end
    end

    return nothing
end

function root_iteration!(root_finder::FixedPoint, dy::Matrix{T}, i::Int64,
                         res::Vector{T}, args...) where T <: AbstractFloat

    @.. dy[:,i] -= res

    return nothing
end

function root_iteration!(root_finder::Newton, dy::Matrix{T}, i::Int64, res::Vector{T},
                         J::Union{Matrix{T}, SparseMatrixCSC{T,Int64}}
                        ) where T <: AbstractFloat

    linear_cache = root_finder.linear_cache
    evaluations = root_finder.evaluations
    time_subroutine = root_finder.time_subroutine
    runtime = root_finder.runtime

    # pass Jacobian and residual error to linear cache
    linear_cache.A = J
    linear_cache.b = res

    if time_subroutine && evaluations[1] % SAMPLE_INTERVAL == 0
        stats = @timed solve!(linear_cache)
        runtime[1] += SAMPLE_INTERVAL*stats.time
    else
        solve!(linear_cache)
    end
    evaluations[1] += 1

    @.. dy[:,i] -= linear_cache.u

    return nothing
end