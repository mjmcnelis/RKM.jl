
# TODO: not really working that well right now, debug later
function evolve_one_time_step!(method::RungeKutta, adaptive::CentralDiff,
                               t::Vector{T}, dt::Vector{T},
                               config::RKMConfig) where T <: AbstractFloat

    update_cache = config.update_cache
    ode_wrap_y! = config.ode_wrap_y!

    iteration = method.iteration
    order = method.order

    epsilon = adaptive.epsilon
    p_norm = adaptive.p_norm
    limiter = adaptive.limiter

    low = limiter.low
    high = limiter.high
    dt_min = limiter.dt_min
    dt_max = limiter.dt_max

    y = update_cache.y
    y_tmp = update_cache.y_tmp
    f_tmp = update_cache.f_tmp
    dy = update_cache.dy

    # note: use y1 to set y_prev
    y_prev = update_cache.y1

    ode_wrap_y!(f_tmp, t[1], y)                             # evaluate first stage at (t,y)

    # TEMP: fixed time step for first update until estimate first time step w/ doubling
    if ode_wrap_y!.FE[1] > 2
        # for high did 1/order, but epsilon 2/x?
        # don't remember my reasoning for that
        # TODO: don't I already rescale these parameters now?
        high    ^= 1.0 / order[1]                       # rescale high based on order
        epsilon ^= 2.0 / (1.0 + order[1])

        # TODO: check if this is allocating (don't set new time step)
        # note: seems like it's allocating even for fixed time step

        # @.. y_tmp = y + dt[1]*f_tmp                   # compute y_star (stored in y_tmp)

        # approximate C w/ central differences (stored in y_tmp)
        # @.. y_tmp = y_tmp - 2.0*y + y_prev
        @.. y_tmp = dt[1]*f_tmp - y + y_prev
        @.. y_tmp *= 2.0/dt[1]^2

        # dy_dt!(y_tmp, t[1], y_prev)
        # @.. y_tmp = f_tmp - y_tmp
        # @.. y_tmp *= 2.0/dt[1]
        # C = 2/dt^2 * (dt*f - y + y_prev)
        # C = 2/dt * (f - (y - y_prev)/dt ) ~ 2/dt * (f_n - f_n-1)

        C_norm = norm(y_tmp, p_norm)
        y_norm = norm(y, p_norm)
        f_norm = norm(f_tmp, p_norm)

        # TODO: solve more complicated algebraic equation from VAH paper
        if iszero(C_norm)                               # compute new time step in dt[2]
            dt[2] = high*dt[1]
        else
            if C_norm*y_norm > 2.0*epsilon*f_norm^2
                dt[2] = sqrt(2.0*epsilon*y_norm/C_norm)
            else
                dt[2] = 2.0*epsilon*f_norm/C_norm
            end
        end
        dt[2] = min(high*dt[1], max(low*dt[1], dt[2]))  # control growth rate
        dt[2] = min(dt_max, max(dt_min, dt[2]))         # impose min/max bounds
    end
    dt[1] = dt[2]                                       # store current time step
    @.. dy[:,1] = dt[1] * f_tmp
    runge_kutta_step!(method, iteration, t, dt, config)

    @.. y_prev = y                                      # store current solution

    return nothing
end