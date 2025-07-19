
function evolve_one_time_step!(method::RungeKutta, adaptive::Fixed,
                               t::Vector{T}, dt::Vector{T},
                               config::RKMConfig) where T <: AbstractFloat

    ode_wrap_y! = config.ode_wrap_y!
    update_cache = config.update_cache
    sensitivity = config.sensitivity
    interpolator = config.interpolator

    iteration = method.iteration
    explicit_stage = method.explicit_stage
    fesal = method.fesal

    dy = update_cache.dy
    y = update_cache.y
    f = update_cache.f
    y_tmp = update_cache.y_tmp
    f_tmp = update_cache.f_tmp

    # evaluate first stage at current time
    if explicit_stage[1]
        @.. dy[:,1] = dt[1] * f
    end

    # for sensitivity
    if explicit_stage[1]
        stage_idx = 1
        @.. y_tmp = y
        explicit_sensitivity_stage!(sensitivity, stage_idx, t[1], dt[1], config, method)
    end

    runge_kutta_step!(method, iteration, t, dt, config)

    # evaluate ODE at next time step and store in f_tmp (skip if method is FESAL)
    if (explicit_stage[1] || interpolator isa CubicHermite) && !fesal
        ode_wrap_y!(f_tmp, t[1] + dt[1], y_tmp)
    end

    return nothing
end