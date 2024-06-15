
function evolve_one_time_step!(method::RungeKutta, adaptive::Doubling,
             controller::Controller, FE::MVector{1,Int64},
             t::Vector{T}, dt::Vector{T},
             ode_wrap!::ODEWrapper, update_cache::RKMCache, linear_cache,
             stage_finder::ImplicitStageFinder) where T <: AbstractFloat

    @unpack epsilon, alpha, delta, p_norm, max_attempts, total_attempts = adaptive
    @unpack explicit_stage = method
    @unpack limiter = controller
    @unpack dt_min, dt_max = limiter
    @unpack y, y_tmp, f, f_tmp, y1, y2, error = update_cache

    order = method.order[1]                             # order of scheme

    if FE[1] == 0
        # always evaluate first stage at initial time (should move outside of function)
        ode_wrap!(f, t[1], y)
        FE[1] += 1
    else
        # get ODE of current time step (should already be stored in f_tmp)
        @.. f = f_tmp
    end

    dt[1] = dt[2]                                       # initialize time step
    rescale = T(1.0)                                    # default time step rescaling

    attempts = 1
    while true                                          # start step doubling routine
        dt[1] = min(dt_max, max(dt_min, dt[1]*rescale)) # increase dt for next attempt

        double_step!(method, t[1], dt[1], ode_wrap!, FE, update_cache,
                     linear_cache, stage_finder)

        @.. error = (y2 - y1) / (2.0^order - 1.0)       # estimate local truncation error
        @.. y2 = y2 + error                             # Richardson extrapolation

        # note: have modified norm function for DoubleFloat
        e_norm = norm(error, p_norm)                # compute norms
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
    ode_wrap!(f_tmp, t[1] + dt[1], y_tmp)
    FE[1] += 1

    return nothing
end

function double_step!(method, t, dt, ode_wrap!, FE, update_cache,
                      linear_cache, stage_finder)

    @unpack explicit_stage, fsal, iteration = method
    @unpack dy, y, y_tmp, f, f_tmp, y1, y2 = update_cache

    # update full time step
    if explicit_stage[1]
        @.. dy[:,1] = dt * f
    end
    runge_kutta_step!(method, iteration, t, dt, ode_wrap!, FE,
                      update_cache, linear_cache, stage_finder)
    @.. y1 = y_tmp

    # update two half time steps
    #   first half step
    if explicit_stage[1]
        @.. dy[:,1] = (dt/2.0) * f
    end
    runge_kutta_step!(method, iteration, t, dt/2, ode_wrap!, FE,
                      update_cache, linear_cache, stage_finder)
    @.. y2 = y_tmp
    #   second half step
    if explicit_stage[1]
        # skip function evaluation if method is FSAL
        if !fsal
            ode_wrap!(f_tmp, t + dt/2.0, y2)
            FE[1] += 1
        end
        @.. dy[:,1] = (dt/2.0) * f_tmp
    end
    # note: swap y, y2 values before/after second runge_kutta_step!
    @.. y_tmp = y2
    @.. y2 = y
    @.. y = y_tmp
    runge_kutta_step!(method, iteration, t+dt/2, dt/2, ode_wrap!, FE,
                      update_cache, linear_cache, stage_finder)
    @.. y = y2
    @.. y2 = y_tmp
    return nothing
end
