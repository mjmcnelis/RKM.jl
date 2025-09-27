
abstract type RootFinderMethod end

struct FixedPoint <: RootFinderMethod
    epsilon::Float64
    p_norm::Float64
    max_iterations::Int64
    iteration_bins::Vector{Int64}
    convergence_failures::MVector{1,Int64}
end

function FixedPoint(; epsilon::Float64 = 1e-8, p_norm::Float64 = 2.0,
                      max_iterations::Int64 = 10)

    @assert p_norm >= 1 "p_norm = $p_norm is not valid"

    iteration_bins = zeros(Int64, max_iterations)
    convergence_failures = MVector{1,Int64}(0)

    return FixedPoint(epsilon, p_norm, max_iterations,
                      iteration_bins, convergence_failures)
end

struct Newton{LC} <: RootFinderMethod where LC <: LinearCache
    linear_cache::LC
    epsilon::Float64
    p_norm::Float64
    max_iterations::Int64
    iteration_bins::Vector{Int64}
    convergence_failures::MVector{1,Int64}
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

    iteration_bins = zeros(Int64, max_iterations)
    convergence_failures = MVector{1,Int64}(0)
    evaluations = MVector{1,Int64}(0)
    time_subroutine = false # dummy default
    runtime = MVector{1,Float64}(0.0)

    return Newton(linear_cache, epsilon, p_norm, max_iterations, iteration_bins,
                  convergence_failures, evaluations, time_subroutine, runtime)
end

function reconstruct_root_finder(root_finder::FixedPoint, res::Vector{T},
                                 J::Union{Matrix{T}, SparseMatrixCSC{T,Int64}},
                                 time_subroutine::Bool) where T <: AbstractFloat

    iteration_bins = root_finder.iteration_bins
    convergence_failures = root_finder.convergence_failures

    iteration_bins .= 0
    convergence_failures[1] = 0

    return root_finder
end

function reconstruct_root_finder(root_finder::Newton, res::Vector{T},
                                 J::Union{Matrix{T}, SparseMatrixCSC{T,Int64}},
                                 time_subroutine::Bool) where T <: AbstractFloat

    iteration_bins = root_finder.iteration_bins
    convergence_failures = root_finder.convergence_failures
    evaluations = root_finder.evaluations
    runtime = root_finder.runtime

    iteration_bins .= 0
    convergence_failures[1] = 0
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