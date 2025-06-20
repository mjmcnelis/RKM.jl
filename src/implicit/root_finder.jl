
abstract type RootFinderMethod end

@kwdef struct Newton{AF} <: RootFinderMethod where AF <: AbstractFactorization
    linear_method::AF = LUFactorization()
    # TODO: reuse adaptive epsilon or 100x smaller?
    epsilon::Float64 = 1e-8
    # TODO: make outer constructor to check p_norm value
    p_norm::Float64 = 2.0
    max_iterations::Int64 = 10
    # TODO: track how many iterations were performed each step
end

@kwdef struct FixedPoint <: RootFinderMethod
    epsilon::Float64 = 1e-8
    p_norm::Float64 = 2.0
    max_iterations::Int64 = 10
    # TODO: track how many iterations were performed each step
end

# TODO: pass dy_stage::SubArray{T} instead
function root_iteration!(root_finder::FixedPoint, dy::Matrix{T}, i::Int64,
                         res::Vector{T}, args...) where T <: AbstractFloat
    @.. dy[:,i] -= res
    return nothing
end

# TODO: pass dy_stage::SubArray{T} instead
function root_iteration!(root_finder::Newton, dy::Matrix{T}, i::Int64, res::Vector{T},
                         J::Union{Matrix{T}, SparseMatrixCSC{T,Int64}},
                         linear_cache) where T <: AbstractFloat
    # pass Jacobian and residual error to linear cache
    linear_cache.A = J
    linear_cache.b = res
    solve!(linear_cache)
    @.. dy[:,i] -= linear_cache.u
    return nothing
end