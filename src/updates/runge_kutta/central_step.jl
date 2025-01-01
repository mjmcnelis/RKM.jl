
# TODO: not really working that well right now, debug later (also out of date)
function evolve_one_time_step!(method::RungeKutta, adaptive::CentralDiff,
             controller::Controller, t::Vector{T}, dt::Vector{T},
             ode_wrap!::ODEWrapperState, update_cache::RKMCache,
             # TODO: may want to do kwargs for different caches used in adaptive methods
             # note: next argument was f but renamed it to y_prev here
             args...) where T <: AbstractFloat

    @unpack iteration = method
    @unpack y, y_tmp, f_tmp, dy = update_cache
    y_prev = update_cache.f

    ode_wrap!(f_tmp, t[1], y)                              # evaluate first stage at (t,y)

    # TEMP: fixed time step for first update until estimate first time step w/ doubling
    if ode_wrap!.FE[1] > 1
        @unpack epsilon, p_norm = adaptive
        @unpack low, high, dt_min, dt_max = controller.limiter

        order = method.order[1]                         # order of scheme

        # for high did 1/order, but epsilon 2/x?
        # don't remember my reasoning for that
        high    ^= 1.0 / order                          # rescale high based on order
        epsilon ^= 2.0 / (1.0 + order)

        # TODO: check if this is allocating (don't set new time step)
        # note: seems like it's allocating even for fixed time step

        # @.. y_tmp = y + dt[1]*f_tmp                     # compute y_star (stored in y_tmp)

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
        if C_norm == 0.0                                # compute new time step in dt[2]
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
    # evaluate first stage iteration w/ new time step (i.e. dt[2])
    @.. dy[:,1] = dt[2] * f_tmp
    runge_kutta_step!(method, iteration, t[1], dt[2], ode_wrap!, update_cache, args...)

    dt[1] = dt[2]                                       # store current time step
    @.. y_prev = y                                      # store current solution
    @.. y = y_tmp                                       # get iteration

    return nothing
end