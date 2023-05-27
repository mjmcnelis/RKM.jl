
abstract type RootMethod end 
struct Newton <: RootMethod end         # TODO: add line search later
struct FixedPoint <: RootMethod end

abstract type JacobianMethod end

@kwdef struct ForwardJacobian{JC} <: JacobianMethod where JC <: JacobianConfig
    cache::JC = JacobianConfig(nothing, [0.0], [0.0])
end

@kwdef struct FiniteJacobian{JC} <: JacobianMethod where JC <: JacobianCache
    cache::JC = JacobianCache([0.0]) 
end

abstract type StageFinder end 

# good enough start (wrap caches later)
@kwdef struct ImplicitStageFinder{JM} <: StageFinder where JM <: JacobianMethod
    root_method::RootMethod = Newton()
    # consider wrapping in Newton struct instead
    jacobian_method::JM     = FiniteJacobian()
    epsilon::Float64        = 1e-8  # TODO: reuse adaptive epsilon or 100x smaller?
    max_iterations::Int64   = 10
    # TODO: make outer constructor to check p_norm value 
    p_norm::Float64         = 2.0
    # add iterations_per_stage
end

function set_jacobian_cache(stage_finder::ImplicitStageFinder, dy_dt!, f, y)
    @unpack jacobian_method = stage_finder
    if jacobian_method isa FiniteJacobian
        cache = JacobianCache(y)
    elseif jacobian_method isa ForwardJacobian 
        cache = JacobianConfig(dy_dt!, f, y)
    end
    @set! stage_finder.jacobian_method.cache = cache
    return stage_finder
end

function evaluate_system_jacobian!(jacobian_method::ForwardJacobian, FE, J, dy_dt!, y, f)
    # TODO: how to reduce allocations here, take it apart or make a wrapper?
    # so in order for jacobian config (i.e. cache) to work, dy_dt! argument
    # has to be the same object stored in jacobian_method.cache
    @unpack cache = jacobian_method 
    jacobian!(J, dy_dt!, f, y, cache)
    FE[1] += ceil(Int64, length(y)/DEFAULT_CHUNK_THRESHOLD)
    return nothing 
end

function evaluate_system_jacobian!(jacobian_method::FiniteJacobian, FE, J, dy_dt!, y, args...)
    @unpack cache = jacobian_method
    finite_difference_jacobian!(J, dy_dt!, y, cache)
    FE[1] += length(y) + 1
    return nothing 
end
