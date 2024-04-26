
function evolve_one_time_step!(method::RungeKutta, adaptive::Embedded,
             controller::Controller, FE::MVector{1,Int64},
             t::Vector{T}, dt::Vector{T},
             ode_wrap!::ODEWrapper, update_cache::RKMCache, linear_cache,
             stage_finder::ImplicitStageFinder) where T <: AbstractFloat

    @unpack epsilon, alpha, delta, p_norm, max_attempts, total_attempts = adaptive
    @unpack iteration, explicit_stage, fsal = method
    @unpack limiter = controller
    @unpack dt_min, dt_max = limiter
    @unpack dy, y, y_tmp, f, f_tmp, y1, y2, error = update_cache

    order_max = maximum(method.order)                   # max/min orders in embedded scheme
    order_min = minimum(method.order)

    # note: comment if benchmark against OrdinaryDiffEq (add boolean?)
    epsilon ^= (order_min / order_max)                  # rescale tolerance parameters
    alpha ^= (order_min / order_max)
    delta ^= (order_min / order_max)
    # TODO: try reconstruct_adaptive (i.e. repower epsilon, make it correct type)
    epsilon = T(epsilon) # tmp for T = Float32 (otherwise tol = Float64 != Float32)
    alpha = T(alpha)
    delta = T(delta)

    # evaluate first stage at (t,y)
    if explicit_stage[1]
        # skip function evaluation if method is FSAL
        if FE[1] > 0 && fsal
            f .= f_tmp
        else
            ode_wrap!(f, t[1], y)
            FE[1] += 1
        end
    end

    dt[1] = dt[2]                                       # initialize time step

    # TODO: float type can change from Float64 to Double64
    # note: is that why benchmark vs OrdinaryDiffEq bad for Double64?
    rescale = 1.0

    attempts = 1
    while true                                          # start embedded routine
        dt[1] = min(dt_max, max(dt_min, dt[1]*rescale)) # increase dt for next attempt

        @.. dy[:,1] = dt[1] * f                         # primary iteration
        runge_kutta_step!(method, iteration, t[1], dt[1], ode_wrap!, FE,
                          update_cache, linear_cache, stage_finder)
        @.. y1 = y_tmp

        embedded_step!(method, y, dy, y_tmp)
        @.. y2 = y_tmp                                  # embedded iteration

        @.. error = y2 - y1                             # local error of embedded pair

        e_norm = norm_tmp(error, p_norm)                    # compute norms
        y_norm = norm_tmp(y1, p_norm)
        # TODO: need to use Δy = y2 since it's secondary method
        #       but should make labeling consistent w/ doubling
        Δy = y2
        @.. Δy = y1 - y
        Δy_norm = norm_tmp(Δy, p_norm)

        # compute tolerance
        tol = max(epsilon*y_norm, alpha, delta*Δy_norm) # compute tolerance

        if !controller.initialized[1]                   # initialize controller
            initialize_controller!(controller, e_norm, tol, dt[1])
        end

        if e_norm == 0.0                                # compute scaling factor for dt
            rescale = limiter.high
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
        attempts <= max_attempts || (println("embedded exceeded $max_attempts attempts");
                                     set_previous_control_vars!(controller, e_norm, tol, dt[1]);
                                     break)
        attempts += 1
    end
    total_attempts[1] += attempts
    controller.initialized[1] = true

    @.. y = y1                                          # get update
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