
function evolve_one_time_step!(method::M, adaptive::Doubling,
                               t::Vector{T}, dt::Vector{T},
                               config::RKMConfig) where {M <: Magnus,
                                                         T <: AbstractFloat}

    ode_wrap_y! = config.ode_wrap_y!
    update_cache = config.update_cache

    epsilon = adaptive.epsilon
    alpha = adaptive.alpha
    delta = adaptive.delta
    p_norm = adaptive.p_norm
    max_attempts = adaptive.max_attempts
    total_attempts = adaptive.total_attempts
    pid = adaptive.pid
    limiter = adaptive.limiter
    initialized_controller = adaptive.initialized_controller

    dt_min = limiter.dt_min
    dt_max = limiter.dt_max

    y = update_cache.y
    y_tmp = update_cache.y_tmp
    f_tmp = update_cache.f_tmp
    y1 = update_cache.y1
    y2 = update_cache.y2
    res = update_cache.res

    order = method.order[1]                             # order of scheme

    dt[1] = dt[2]                                       # initialize time step
    rescale = T(1.0)                                    # default time step rescaling

    attempts = 1
    while true                                          # start step doubling routine
        dt[1] = min(dt_max, max(dt_min, dt[1]*rescale)) # increase dt for next attempt

        double_step!(method, t, dt, config)

        @.. res = (y2 - y1) / (2.0^order - 1.0)         # estimate local truncation error

        @.. y2 = y2 + res

        # note: have modified norm function for DoubleFloat
        e_norm = norm(res, p_norm)                      # compute norms
        y_norm = norm(y2, p_norm)
        Δy = y1
        @.. Δy = y2 - y
        Δy_norm = norm(Δy, p_norm)

        tol = max(epsilon*y_norm, alpha, delta*Δy_norm) # compute tolerance

        if !initialized_controller[1]                   # initialize controller variables
            initialize_controller!(update_cache, e_norm, tol, dt[1])
        end

        if iszero(e_norm)                               # compute scaling factor for dt
            rescale = T(limiter.high)
        else
            rescale = rescale_time_step(pid, update_cache, tol, e_norm)
            rescale = limit_time_step(limiter, rescale)
            # TODO: need to track rescale in the controller
            #       both for accepted and rejected step
        end

        dt[2] = min(dt_max, max(dt_min, dt[1]*rescale)) # projected dt for next iteration

        if e_norm <= tol                                # compare error to tolerance
            set_previous_control_vars!(update_cache, e_norm, tol, dt[1])
            break
        end
        attempts <= max_attempts || (println("step doubling exceeded $max_attempts attempts");
                                     set_previous_control_vars!(update_cache, e_norm, tol, dt[1]);
                                     break)
        attempts += 1
    end
    total_attempts[1] += attempts
    initialized_controller[1] = true

    @.. y_tmp = y2                                      # get iteration

    # evaluate ODE at next time step and store in f_tmp
    # note: if do Richardson extrapolation, then always have
    # to evaluate (explicit) first stage at (t+dt,y+dy)
    ode_wrap_y!(f_tmp, t[1] + dt[1], y_tmp)

    return nothing
end

function double_step!(method::M, t::Vector{T},
                      dt::Vector{T}, config::RKMConfig) where {M <: Magnus,
                                                               T <: AbstractFloat}
    update_cache = config.update_cache

    y = update_cache.y
    y_tmp = update_cache.y_tmp
    y1 = update_cache.y1
    y2 = update_cache.y2

    t_start = t[1]
    dt_start = dt[1]

    # update two half time steps
    dt[1] /= 2.0

    # first half step
    magnus_step!(method, t, dt, config)
    @.. y2 = y_tmp

    # second half step
    # note: swap y, y2 values before/after second magnus_step!
    t[1] += dt[1]
    @.. y_tmp = y2
    @.. y2 = y
    @.. y = y_tmp
    magnus_step!(method, t, dt, config)
    @.. y = y2
    @.. y2 = y_tmp

    # undo changes to t, dt after double step finished
    t[1] = t_start
    dt[1] = dt_start

    # update full time step
    magnus_step!(method, t, dt, config)
    @.. y1 = y_tmp

    return nothing
end
