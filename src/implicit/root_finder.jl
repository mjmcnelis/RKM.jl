
abstract type RootFinderMethod end

@kwdef struct FixedPoint <: RootFinderMethod
    epsilon::Float64 = 1e-8
    p_norm::Float64 = 2.0
    max_iterations::Int64 = 10
    # TODO: track how many iterations were performed each step
end

struct Newton{LC} <: RootFinderMethod where LC <: LinearCache
    linear_cache::LC
    # TODO: reuse adaptive epsilon or 100x smaller?
    epsilon::Float64
    # TODO: make outer constructor to check p_norm value
    p_norm::Float64
    max_iterations::Int64
    # TODO: track how many iterations were performed each step
end

function Newton(; linear_method::AF = LUFactorization(),
                  epsilon::Float64 = 1e-8, p_norm::Float64 = 2.0,
                  max_iterations::Int64 = 10) where AF <: AbstractFactorization

    # configure linear cache (see src/common.jl in LinearSolve.jl)
    J = zeros(0, 0)
    res = zeros(0)
    linear_cache = init(LinearProblem(J, res), linear_method,
                        alias = LinearAliasSpecifier(; alias_A = true, alias_b = true))

    return Newton(linear_cache, epsilon, p_norm, max_iterations)
end

function reconstruct_root_finder(root_finder::FixedPoint, args...)
    return root_finder
end

function reconstruct_root_finder(root_finder::Newton, res::Vector{T},
             J::Union{Matrix{T}, SparseMatrixCSC{T,Int64}}) where T <: AbstractFloat

    @unpack linear_cache = root_finder
    linear_method = linear_cache.alg

    linear_cache = init(LinearProblem(J, res), linear_method;
                        alias = LinearAliasSpecifier(; alias_A = true, alias_b = true),)
    @set! root_finder.linear_cache = linear_cache

    return root_finder
end

function root_iteration!(root_finder::FixedPoint, dy_stage::SubArray{T},
                         res::Vector{T}, args...) where T <: AbstractFloat
    @.. dy_stage -= res
    return nothing
end

function root_iteration!(root_finder::Newton, dy_stage::SubArray{T}, res::Vector{T},
                         J::Union{Matrix{T}, SparseMatrixCSC{T,Int64}}
                        ) where T <: AbstractFloat

    @unpack linear_cache = root_finder
    # pass Jacobian and residual error to linear cache
    linear_cache.A = J
    linear_cache.b = res
    solve!(linear_cache)
    @.. dy_stage -= linear_cache.u

    return nothing
end