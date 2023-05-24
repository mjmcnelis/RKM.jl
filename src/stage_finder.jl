
abstract type RootMethod end 
struct Newton <: RootMethod end         # TODO: add line search later
struct FixedPoint <: RootMethod end

abstract type JacobianMethod end
struct ForwardJacobian <: JacobianMethod end

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
    # add iterations_per_stage, p_norm
end

function set_jacobian_cache(stage_finder::ImplicitStageFinder, y)
    if stage_finder.jacobian_method isa FiniteJacobian
        @set! stage_finder.jacobian_method.cache = JacobianCache(y)
    end
    return stage_finder
end

function evaluate_system_jacobian!(jacobian_method::ForwardJacobian, J, dy_dt!, f, y,
                                   jacobian_config)
    # TODO: how to reduce allocations here, take it apart or make a wrapper?
    # so in order for jacobian config to work, dy_dt_wrap! argument
    # has to be the same object stored in jacobian_config
    jacobian!(J, dy_dt!, f, y, jacobian_config)
    return nothing 
end

function evaluate_system_jacobian!(jacobian_method::FiniteJacobian, J, dy_dt!, f, y, 
                                   jacobian_config)
    @unpack cache = jacobian_method
    finite_difference_jacobian!(J, dy_dt!, y, cache)
    return nothing 
end
