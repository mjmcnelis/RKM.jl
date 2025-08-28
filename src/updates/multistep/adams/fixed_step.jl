
function evolve_one_time_step!(method::Adams, adaptive::Fixed,
                               t::Vector{T}, dt::Vector{T},
                               config::RKMConfig) where T <: AbstractFloat

    ode_wrap_y! = config.ode_wrap_y!
    update_cache = config.update_cache

    stages = method.stages
    iteration = method.iteration
    order = method.order
    start_counter = method.start_counter

    dy = update_cache.dy
    dy_LM = update_cache.dy_LM
    y = update_cache.y
    y_tmp = update_cache.y_tmp
    f = update_cache.f

    # evaluate ODE at current state/time (i.e E)
    ode_wrap_y!(f, t[1], y)

    @.. dy[:,1] = dt[1] * f
    @.. dy_LM[:,1] = dt[1] * f

    # initialize previous stages (trivial)
    # TODO: check if this works (moved FE to ode_wrap!)
    # if ode_wrap!.evaluations[1] == 1
    #     for j in 2:stages
    #         @.. dy_LM[:,j] = dt[1] * f
    #     end
    # end

    if start_counter[1] < order-1
        start_method = method.start_method
        start_iteration = start_method.iteration

        runge_kutta_step!(start_method, start_iteration, t, dt, config)

        start_counter[1] += 1
    else
        adams_step!(method, iteration, t, dt, config)
    end

    # shift stages (except first one)
    for j in 1:stages-1
        @.. dy_LM[:,stages+1-j] = dy_LM[:,stages-j]
    end
    @.. y = y_tmp
    return nothing
end

function evolve_one_time_step!(method::DifferentiationFormula, adaptive::Fixed,
             t::Vector{T}, dt::Vector{T}, ode_wrap!::ODEWrapperState,
             update_cache::UpdateCache, linear_cache,
             state_jacobian::JacobianMethod) where T <: AbstractFloat

    stages = method.stages
    iteration = method.iteration

    dy = update_cache.dy
    y = update_cache.y
    y_tmp = update_cache.y_tmp

    # initialize previous states (trivial)
    # TODO: check if this works (moved FE to ode_wrap!)
    if ode_wrap!.evaluations[1] == 1
        for j in 2:stages
            @.. dy[:,j] = y
        end
    end

    differentiation_formula_step!(method, iteration, t[1], dt[1], ode_wrap!,
                                  update_cache, linear_cache, state_jacobian)

    # shift states
    for j in 1:stages-2
        # @show stages+1-j, stages-j
        @.. dy[:,stages+1-j] = dy[:,stages-j]
    end
    @.. dy[:,2] = y

    # get update
    @.. y = y_tmp
    return nothing
end