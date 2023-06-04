
abstract type RootMethod end 
# struct Newton <: RootMethod end
@kwdef struct Newton <: RootMethod 
    switch_broyden::Bool = true 
end         # TODO: add line search later

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
@kwdef struct ImplicitStageFinder{JM} <: StageFinder where JM <: JacobianMethod
    root_method::RootMethod = Newton()
    jacobian_method::JM     = FiniteJacobian()
     # TODO: reuse adaptive epsilon or 100x smaller?
    epsilon::Float64        = 1e-8 
    max_iterations::Int64   = 10
    # TODO: make outer constructor to check p_norm value 
    p_norm::Float64         = 2.0
    # add iterations_per_stage
    iterations::Vector{UInt16} = UInt16[]
end

function set_jacobian_cache(stage_finder::ImplicitStageFinder, dy_dt!, f, y)
    # TODO: only set cache is use Newton method
    @unpack jacobian_method = stage_finder
    if jacobian_method isa FiniteJacobian
        cache = JacobianCache(y)
    elseif jacobian_method isa ForwardJacobian 
        cache = JacobianConfig(dy_dt!, f, y)
    end
    @set! stage_finder.jacobian_method.cache = cache
    return stage_finder
end

function evaluate_system_jacobian!(jacobian_method::ForwardJacobian, 
                                   root_method::Newton, n,
                                   FE, J, dy_dt!, y, dt, A, f)
    # TODO: how to reduce allocations here, take it apart or make a wrapper?
    # so in order for jacobian config (i.e. cache) to work, dy_dt! argument
    # has to be the same object stored in jacobian_method.cache
    @unpack cache, evaluations = jacobian_method 
    @unpack switch_broyden = root_method 

    if n > 1 && switch_broyden
        # TODO: wrap function 
        dx2 = dot(linear_cache.u, linear_cache.u)
        for i in 1:2 
            for k = 1:2 
                error_prev[i] -= J[i,k] * linear_cache.u[k]
            end 
        end
        for i in 1:2
            for j in 1:2
                J[i,j] -= (error[i] - error_prev[i]) * linear_cache.u[j] / dx2
            end
        end
        @show J 
        println("")
    else
        jacobian!(J, dy_dt!, f, y, cache)
        J .*= -A*dt                                         # J <- I - A.dt.J
        for i in diagind(J)
            J[i] += 1.0
        end
    end
    FE[1] += ceil(Int64, length(y)/DEFAULT_CHUNK_THRESHOLD)
    evaluations[1] += 1
    return nothing 
end

function evaluate_system_jacobian!(jacobian_method::FiniteJacobian, 
                                   root_method::Newton, n,
                                   FE, J, dy_dt!, y, dt, A, args...)
    @unpack cache, evaluations = jacobian_method
    @unpack switch_broyden = root_method

    if n > 1 && switch_broyden
        # TODO: wrap function 
        dx2 = dot(linear_cache.u, linear_cache.u)
        for i in 1:2 
            for k = 1:2 
                error_prev[i] -= J[i,k] * linear_cache.u[k]
            end 
        end
        for i in 1:2
            for j in 1:2
                J[i,j] -= (error[i] - error_prev[i]) * linear_cache.u[j] / dx2
            end
        end
        @show J 
        println("")
    else
        finite_difference_jacobian!(J, dy_dt!, y, cache)
        J .*= -A*dt                                         # J <- I - A.dt.J
        for i in diagind(J)
            J[i] += 1.0
        end
    end

    FE[1] += length(y) + 1
    evaluations[1] += 1
    return nothing 
end
