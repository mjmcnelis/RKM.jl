
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
    # reset jacobian evaluations
    stage_finder.jacobian_method.evaluations .= 0

    # TODO: only set cache is use Newton method
    @unpack jacobian_method = stage_finder
    if jacobian_method isa FiniteJacobian
        @unpack sparsity = jacobian_method
        if all(size(sparsity) .== length(y))
            colorvec = matrix_colors(sparsity)
            cache = JacobianCache(y; colorvec, sparsity)
        else
            cache = JacobianCache(y)
        end
    elseif jacobian_method isa ForwardJacobian
        cache = JacobianConfig(ode_wrap!, f, y)
    elseif jacobian_method isa ForwardColorJacobian
        @unpack sparsity = jacobian_method
        if all(size(sparsity) .== length(y))
            colorvec = matrix_colors(sparsity)
            cache = ForwardColorJacCache(ode_wrap!, y; colorvec, sparsity)
        else
            cache = ForwardColorJacCache(ode_wrap!, y)
        end
    end
    @set! stage_finder.jacobian_method.cache = cache
    return stage_finder
end
