
abstract type RootMethod end 
struct Newton <: RootMethod end         # TODO: add line search later
struct FixedPoint <: RootMethod end

abstract type JacobianMethod end
struct ForwardJacobian <: JacobianMethod end 
struct FiniteJacobian <: JacobianMethod end

abstract type StageFinder end 

# good enough start (wrap caches later)
@kwdef struct ImplicitStageFinder <: StageFinder
    root_method::RootMethod         = Newton()
    # consider wrapping in Newton struct instead
    jacobian_method::JacobianMethod = FiniteJacobian()
    epsilon::Float64                = 1e-8  # TODO: reuse adaptive epsilon or 100x smaller?
    max_iterations::Int64           = 10
    # add iterations_per_stage, p_norm
end

function evaluate_system_jacobian!(::ForwardJacobian, J, dy_dt!, f, y,
                                   jacobian_config, finitediff_cache)
    # TODO: how to reduce allocations here, take it apart or make a wrapper?
    # so in order for jacobian config to work, dy_dt_wrap! argument
    # has to be the same object stored in jacobian_config
    jacobian!(J, dy_dt!, f, y, jacobian_config)
    return nothing 
end

function evaluate_system_jacobian!(::FiniteJacobian, J, dy_dt!, f, 
                                   y, jacobian_config, finitediff_cache)
    finite_difference_jacobian!(J, dy_dt!, y, finitediff_cache)
    return nothing 
end