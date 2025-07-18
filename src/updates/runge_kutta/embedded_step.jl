
function evolve_one_time_step!(method::RungeKutta, adaptive::Embedded,
                               t::Vector{T}, dt::Vector{T},
                               config::RKMConfig) where T <: AbstractFloat

    ode_wrap_y! = config.ode_wrap_y!
    update_cache = config.update_cache
    interpolator = config.interpolator

    iteration = method.iteration
    explicit_stage = method.explicit_stage
    fesal = method.fesal

    epsilon = adaptive.epsilon
    alpha = adaptive.alpha
    delta = adaptive.delta
    p_norm = adaptive.p_norm
    max_attempts = adaptive.max_attempts
    total_attempts = adaptive.total_attempts
    limiter = adaptive.limiter
    initialized_controller = adaptive.initialized_controller

    dt_min = limiter.dt_min
    dt_max = limiter.dt_max

    dy = update_cache.dy
    y = update_cache.y
    y_tmp = update_cache.y_tmp
    f = update_cache.f
    f_tmp = update_cache.f_tmp
    y1 = update_cache.y1
    y2 = update_cache.y2
    res = update_cache.res

    dt[1] = dt[2]                                       # initialize time step
    rescale = T(1.0)                                    # default time step rescaling

    attempts = 1
    while true                                          # start embedded routine
        dt[1] = min(dt_max, max(dt_min, dt[1]*rescale)) # increase dt for next attempt

        if explicit_stage[1]
            @.. dy[:,1] = dt[1] * f
        end
        runge_kutta_step!(method, iteration, t, dt, config)
        @.. y1 = y_tmp                                  # primary iteration

        embedded_step!(method, y, dy, y_tmp)
        @.. y2 = y_tmp                                  # embedded iteration

        @.. res = y2 - y1                               # local error of embedded pair

        e_norm = norm(res, p_norm)                      # compute norms
        y_norm = norm(y1, p_norm)
        # TODO: need to use Δy = y2 since it's secondary method
        #       but should make labeling consistent w/ doubling
        Δy = y2
        @.. Δy = y1 - y
        Δy_norm = norm(Δy, p_norm)

        # compute tolerance
        tol = max(epsilon*y_norm, alpha, delta*Δy_norm) # compute tolerance

        if !initialized_controller[1]                   # initialize controller variables
            initialize_controller!(update_cache, e_norm, tol, dt[1])
        end

        if e_norm == 0.0                                # compute scaling factor for dt
            rescale = T(limiter.high)
        else
            # TODO: pass pid instead of adaptive
            rescale = rescale_time_step(adaptive, update_cache, tol, e_norm)
            rescale = limit_time_step(limiter, rescale)
            # TODO: need to track rescale in the controller
            #       both for accepted and rejected step
        end

        dt[2] = min(dt_max, max(dt_min, dt[1]*rescale)) # projected dt for next iteration

        if e_norm <= tol                                # compare error to tolerance
            set_previous_control_vars!(update_cache, e_norm, tol, dt[1])
            break
        end
        attempts <= max_attempts || (println("embedded exceeded $max_attempts attempts");
                                     set_previous_control_vars!(update_cache, e_norm, tol, dt[1]);
                                     break)
        attempts += 1
    end
    total_attempts[1] += attempts
    initialized_controller[1] = true

    @.. y_tmp = y1                                      # get iteration

    # evaluate ODE at next time step and store in f_tmp (skip if method is FESAL)
    if (explicit_stage[1] || interpolator isa CubicHermite) && !fesal
        ode_wrap_y!(f_tmp, t[1] + dt[1], y_tmp)
    end

    return nothing
end

@muladd function embedded_step!(method, y, dy, y_tmp)

    stages = method.stages
    b_hat = method.b_hat

    @.. y_tmp = y                                       # evaluate iteration
    for j in 1:stages
        @.. y_tmp = y_tmp + b_hat[j]*dy[:,j]
    end

    return nothing
end