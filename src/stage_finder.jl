
abstract type RootMethod end
struct Newton <: RootMethod end         # TODO: add line search later
struct FixedPoint <: RootMethod end

abstract type JacobianMethod end

@kwdef struct ForwardJacobian{JC} <: JacobianMethod where JC <: JacobianConfig
    cache::JC = JacobianConfig(nothing, [0.0], [0.0])
    evaluations::MVector{1, Int64} = MVector{1,Int64}(0)
end

@kwdef struct FiniteJacobian{JC} <: JacobianMethod where JC <: JacobianCache
    cache::JC = JacobianCache([0.0])
    evaluations::MVector{1, Int64} = MVector{1,Int64}(0)
end

abstract type StageFinder end

# good enough start (wrap caches later)
@kwdef struct ImplicitStageFinder{JM, AF} <: StageFinder where {JM <: JacobianMethod,
                                                                AF <: AbstractFactorization
                                                               }
    root_method::RootMethod = Newton()
    jacobian_method::JM     = FiniteJacobian()
    linear_method::AF       = LUFactorization()
     # TODO: reuse adaptive epsilon or 100x smaller?
    epsilon::Float64        = 1e-8
    max_iterations::Int64   = 10
    # TODO: make outer constructor to check p_norm value
    p_norm::Float64         = 2.0
    # add iterations_per_stage
end

function set_jacobian_cache(stage_finder::ImplicitStageFinder, ode_wrap!, f, y)
    # TODO: only set cache is use Newton method
    @unpack jacobian_method = stage_finder
    if jacobian_method isa FiniteJacobian
        cache = JacobianCache(y)
    elseif jacobian_method isa ForwardJacobian
        cache = JacobianConfig(ode_wrap!, f, y)
    end
    @set! stage_finder.jacobian_method.cache = cache
    return stage_finder
end

function evaluate_system_jacobian!(jacobian_method::ForwardJacobian, FE, J, ode_wrap!, y, f)
    @unpack cache, evaluations = jacobian_method
    jacobian!(J, ode_wrap!, f, y, cache)
    FE[1] += ceil(Int64, length(y)/DEFAULT_CHUNK_THRESHOLD)
    evaluations[1] += 1
    return nothing
end

function evaluate_system_jacobian!(jacobian_method::FiniteJacobian, FE, J, ode_wrap!, y, args...)
    @unpack cache, evaluations = jacobian_method
    finite_difference_jacobian!(J, ode_wrap!, y, cache)
    FE[1] += length(y) + 1
    evaluations[1] += 1
    return nothing
end
