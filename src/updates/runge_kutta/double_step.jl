
function evolve_one_time_step!(method::RungeKutta, adaptive::Doubling,
             controller::Controller, t::Vector{T}, dt::Vector{T},
             ode_wrap_y!::ODEWrapperState, update_cache::RKMCache, linear_cache,
             stage_finder::ImplicitStageFinder,
             # note: sensitivity not implemented for double step yet
             sensitivity::SensitivityMethod, ode_wrap_p!::ODEWrapperParam,
             interpolator::Interpolator) where T <: AbstractFloat

    @unpack epsilon, alpha, delta, p_norm, max_attempts, total_attempts = adaptive
    @unpack explicit_stage = method
    @unpack limiter = controller
    @unpack dt_min, dt_max = limiter
    @unpack y, y_tmp, f, f_tmp, y1, y2, res = update_cache

    order = method.order[1]                             # order of scheme

    dt[1] = dt[2]                                       # initialize time step
    rescale = T(1.0)                                    # default time step rescaling

    attempts = 1
    while true                                          # start step doubling routine
        dt[1] = min(dt_max, max(dt_min, dt[1]*rescale)) # increase dt for next attempt

        double_step!(method, t, dt, ode_wrap_y!, update_cache, linear_cache,
                     stage_finder, sensitivity, ode_wrap_p!)

        @.. res = (y2 - y1) / (2.0^order - 1.0)     # estimate local truncation error
        @.. y2 = y2 + res                           # Richardson extrapolation

        # note: have modified norm function for DoubleFloat
        e_norm = norm(res, p_norm)                  # compute norms
        y_norm = norm(y2, p_norm)
        Δy = y1
        @.. Δy = y2 - y
        Δy_norm = norm(Δy, p_norm)

        tol = max(epsilon*y_norm, alpha, delta*Δy_norm) # compute tolerance

        if !controller.initialized[1]                   # initialize controller
            initialize_controller!(controller, e_norm, tol, dt[1])
        end

        if e_norm == 0.0                                # compute scaling factor for dt
            rescale = T(limiter.high)
        else
            rescale = rescale_time_step(controller, tol, e_norm)
            rescale = limit_time_step(limiter, rescale)
            # TODO: need to track rescale in the controller
            #       both for accepted and rejected step
        end

        dt[2] = min(dt_max, max(dt_min, dt[1]*rescale)) # projected dt for next iteration

        if e_norm <= tol                                # compare error to tolerance
            set_previous_control_vars!(controller, e_norm, tol, dt[1])
            break
        end
        attempts <= max_attempts || (println("step doubling exceeded $max_attempts attempts");
                                     set_previous_control_vars!(controller, e_norm, tol, dt[1]);
                                     break)
        attempts += 1
    end
    total_attempts[1] += attempts
    controller.initialized[1] = true

    @.. y_tmp = y2                                      # get iteration

    # evaluate ODE at next time step and store in f_tmp
    # note: if do Richardson extrapolation, then always have
    # to evaluate (explicit) first stage at (t+dt,y+dy)
    ode_wrap_y!(f_tmp, t[1] + dt[1], y_tmp)

    return nothing
end

function double_step!(method, t, dt, ode_wrap_y!, update_cache, linear_cache,
                      stage_finder, sensitivity, ode_wrap_p!)

    @unpack explicit_stage, fesal, iteration = method
    @unpack dy, y, y_tmp, f, f_tmp, y1, y2 = update_cache

    t_start = t[1]
    dt_start = dt[1]

    # update full time step
    if explicit_stage[1]
        @.. dy[:,1] = dt[1] * f
    end
    runge_kutta_step!(method, iteration, t, dt, ode_wrap_y!, update_cache,
                      linear_cache, stage_finder, sensitivity, ode_wrap_p!)
    @.. y1 = y_tmp

    # update two half time steps
    dt[1] /= 2.0
    #   first half step
    if explicit_stage[1]
        @.. dy[:,1] = dt[1] * f
    end

    runge_kutta_step!(method, iteration, t, dt, ode_wrap_y!, update_cache,
                      linear_cache, stage_finder, sensitivity, ode_wrap_p!)
    @.. y2 = y_tmp
    #   second half step
    t[1] += dt[1]
    if explicit_stage[1]
        # skip function evaluation if method is FESAL
        if !fesal
            ode_wrap_y!(f_tmp, t[1], y2)
        end
        @.. dy[:,1] = dt[1] * f_tmp
    end
    # note: swap y, y2 values before/after second runge_kutta_step!
    @.. y_tmp = y2
    @.. y2 = y
    @.. y = y_tmp
    runge_kutta_step!(method, iteration, t, dt, ode_wrap_y!, update_cache,
                      linear_cache, stage_finder, sensitivity, ode_wrap_p!)
    @.. y = y2
    @.. y2 = y_tmp

    # undo changes to t, dt after double step finished
    t[1] = t_start
    dt[1] = dt_start
    return nothing
end
