
function evolve_one_time_step!(method::RungeKutta, iteration::Iteration,
             adaptive::Doubling, controller::Controller, FE::MVector{1,Int64},  
             y::Vector{T}, t::Union{Vector{T}, MVector{1,T}}, 
             dt::Union{Vector{T}, MVector{2,T}}, dy_dt!::F, dy::Matrix{T}, 
             y_tmp::Vector{T}, f_tmp::Vector{T}, f::Vector{T}, y1::Vector{T}, 
             y2::Vector{T}, error::Vector{T}, 
             J::MatrixMMatrix, linear_cache, dy_dt_wrap!::ODEWrapper, 
             stage_finder::ImplicitStageFinder) where {T <: AbstractFloat, F}

    @unpack epsilon, low, high, safety, p_norm, dt_min, dt_max, max_attempts = adaptive

    order = method.order[1]                             # order of scheme

    # note: comment if want to compare to OrdinaryDiffEq
    high ^= order / (1.0 + order)                       # rescale high based on order
    epsilon ^= order / (1.0 + order)

    if iteration isa Explicit
        dy_dt!(f, t[1], y)                              # evaluate first stage at (t,y)
        FE[1] += 1
    end

    dt[1] = dt[2]                                       # initialize time step
    rescale = 1.0

    attempts = 1
    while true                                          # start step doubling routine
        dt[1] = min(dt_max, max(dt_min, dt[1]*rescale)) # increase dt for next attempt

        double_step!(method, iteration, y, t[1], dt[1], dy_dt!, dy, y_tmp, f_tmp, 
                     f, y1, y2, FE, J, linear_cache, dy_dt_wrap!, stage_finder)

        @.. error = (y2 - y1) / (2.0^order - 1.0)       # estimate local truncation error
        @.. y2 = y2 + error                             # Richardson extrapolation

        e_norm = norm(error, p_norm)                    # compute norms
        y_norm = norm(y2, p_norm)
        Δy = y1
        @.. Δy = y2 - y
        Δy_norm = norm(Δy, p_norm)

        tol = epsilon * max(y_norm, Δy_norm)            # compute tolerance

        if !controller.initialized[1]                   # initialize controller
            initialize_controller!(controller, e_norm, dt[1])
        end

        if e_norm == 0.0                                # compute scaling factor for dt
            rescale = high
        else
            rescale = rescale_time_step(controller, tol, e_norm, order)
            rescale = min(high, max(low, safety*rescale))
        end

        dt[2] = min(dt_max, max(dt_min, dt[1]*rescale)) # projected dt for next iteration

        if e_norm <= tol                                # compare error to tolerance
            set_previous_control_vars!(controller, e_norm, dt[1])
            break 
        end
        attempts <= max_attempts || (@warn "step doubling exceeded $max_attempts attempts"; break)
        attempts += 1
    end
    controller.initialized[1] = true

    @.. y = y2                                          # get iteration
    return nothing
end

function double_step!(method, iteration, y, t, dt, dy_dt!, dy, y_tmp, f_tmp,
                      f, y1, y2, FE, J, linear_cache, dy_dt_wrap!, stage_finder)

    # note: dy = dt * f lines only matter for explicit or ESDIRK methods
    #       in implicit version, those lines and dy_dt! weren't needed
    # iterate full time step
    @.. dy[:,1] = dt * f
    runge_kutta_step!(method, iteration, y, t, dt, dy_dt!, dy, y_tmp, 
                      f_tmp, FE, J, linear_cache, dy_dt_wrap!, stage_finder)
    @.. y1 = y_tmp

    # iterate two half time steps
    @.. dy[:,1] = (dt/2.0) * f
    runge_kutta_step!(method, iteration, y, t, dt/2., dy_dt!, dy, y_tmp, 
                      f_tmp, FE, J, linear_cache, dy_dt_wrap!, stage_finder)
    @.. y2 = y_tmp
    # note: for some reason, dy_dt! step reduces allocations in implicit
    # TODO: always evaluate stage at t + c[1]*dt (if explicit)
    dy_dt!(f_tmp, t + dt/2.0, y2)#; skip = false)
    FE[1] += 1
    @.. dy[:,1] = (dt/2.0) * f_tmp
    runge_kutta_step!(method, iteration, y2, t+dt/2.0, dt/2.0, dy_dt!, dy, y_tmp,
                      f_tmp, FE, J, linear_cache, dy_dt_wrap!, stage_finder)
    @.. y2 = y_tmp
    return nothing
end
