
function evolve_one_time_step!(method::RungeKutta, adaptive::Embedded,
             t::Vector{T}, dt::Vector{T}, ode_wrap_y!::ODEWrapperState,
             update_cache::RKMCache, state_jacobian::JacobianMethod,
             root_finder::RootFinderMethod, eigenmax::EigenMaxMethod,
             # note: sensitivity not implemented for embedded step yet
             sensitivity::SensitivityMethod, ode_wrap_p!::ODEWrapperParam,
             interpolator::Interpolator) where T <: AbstractFloat

    @unpack epsilon, alpha, delta, p_norm, max_attempts,
            total_attempts, limiter, initialized_controller = adaptive
    @unpack iteration, explicit_stage, fesal = method
    @unpack dt_min, dt_max = limiter
    @unpack dy, y, y_tmp, f, f_tmp, y1, y2, res = update_cache

    dt[1] = dt[2]                                       # initialize time step
    rescale = T(1.0)                                    # default time step rescaling

    attempts = 1
    while true                                          # start embedded routine
        dt[1] = min(dt_max, max(dt_min, dt[1]*rescale)) # increase dt for next attempt

        if explicit_stage[1]
            @.. dy[:,1] = dt[1] * f
        end
        runge_kutta_step!(method, iteration, t, dt, ode_wrap_y!, update_cache, state_jacobian,
                          root_finder, eigenmax, sensitivity, ode_wrap_p!)
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
    @unpack stages, b_hat = method
    @.. y_tmp = y                                       # evaluate iteration
    for j in 1:stages
        dy_stage = view(dy,:,j)
        @.. y_tmp = y_tmp + b_hat[j]*dy_stage
    end
    return nothing
end