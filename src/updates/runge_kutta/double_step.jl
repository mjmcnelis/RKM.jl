
function evolve_one_time_step!(method::RungeKutta, iteration::Iteration,
             adaptive::Doubling, controller::Controller, FE::MVector{1,Int64},
             y::Vector{T}, t::Union{Vector{T}, MVector{1,T}},
             dt::Union{Vector{T}, MVector{2,T}}, ode_wrap!::ODEWrapper, dy::Matrix{T},
             y_tmp::Vector{T}, f_tmp::Vector{T}, f::Vector{T}, y1::Vector{T},
             y2::Vector{T}, error::Vector{T},
             J::MatrixMMatrix, linear_cache,
             stage_finder::ImplicitStageFinder) where T <: AbstractFloat

    @unpack epsilon, low, high, safety, p_norm, dt_min, dt_max, max_attempts = adaptive
    @unpack explicit_stage, fsal = method

    order = method.order[1]                             # order of scheme

    # note: comment if want to compare to OrdinaryDiffEq
    high ^= order / (1.0 + order)                       # rescale high based on order
    epsilon ^= order / (1.0 + order)

    # if do Richardson extrapolation, then always have
    # to evaluate (explicit) first stage at (t,y)
    if explicit_stage[1]
        ode_wrap!(f, t[1], y)
        FE[1] += 1
    end

    dt[1] = dt[2]                                       # initialize time step
    rescale = 1.0

    attempts = 1
    while true                                          # start step doubling routine
        dt[1] = min(dt_max, max(dt_min, dt[1]*rescale)) # increase dt for next attempt

        double_step!(method, iteration, y, t[1], dt[1], ode_wrap!, dy, y_tmp, f_tmp,
                     f, y1, y2, FE, error, J, linear_cache, stage_finder)

        @.. error = (y2 - y1) / (2.0^order - 1.0)       # estimate local truncation error
        @.. y2 = y2 + error                             # Richardson extrapolation

        e_norm = norm(error, p_norm)                    # compute norms
        y_norm = norm(y2, p_norm)
        Δy = y1
        @.. Δy = y2 - y
        Δy_norm = norm(Δy, p_norm)

        tol = epsilon * max(y_norm, Δy_norm)            # compute tolerance

        if !controller.initialized[1]                   # initialize controller
            initialize_controller!(controller, e_norm, tol, dt[1])
        end

        if e_norm == 0.0                                # compute scaling factor for dt
            @warn "Error estimate is zero"
            rescale = high
        else
            rescale = rescale_time_step(controller, tol, e_norm)
            rescale = min(high, max(low, safety*rescale))
        end

        dt[2] = min(dt_max, max(dt_min, dt[1]*rescale)) # projected dt for next iteration

        if e_norm <= tol                                # compare error to tolerance
            set_previous_control_vars!(controller, e_norm, tol, dt[1])
            break
        end
        attempts <= max_attempts || (@warn "step doubling exceeded $max_attempts attempts";
                                     set_previous_control_vars!(controller, e_norm, tol, dt[1]);
                                     break)
        attempts += 1
    end
    controller.initialized[1] = true

    @.. y = y2                                          # get iteration
    return nothing
end

function double_step!(method, iteration, y, t, dt, ode_wrap!, dy, y_tmp, f_tmp,
                      f, y1, y2, FE, error, J, linear_cache, stage_finder)

    @unpack explicit_stage, fsal = method

    # update full time step
    if explicit_stage[1]
        @.. dy[:,1] = dt * f
    end
    runge_kutta_step!(method, iteration, y, t, dt, ode_wrap!, dy, y_tmp,
                      f_tmp, FE, error, J, linear_cache, stage_finder)
    @.. y1 = y_tmp

    # update two half time steps
    #   first half step
    if explicit_stage[1]
        @.. dy[:,1] = (dt/2.0) * f
    end
    runge_kutta_step!(method, iteration, y, t, dt/2., ode_wrap!, dy, y_tmp,
                      f_tmp, FE, error, J, linear_cache, stage_finder)
    @.. y2 = y_tmp
    #   second half step
    if explicit_stage[1]
        # skip function evaluation if method is FSAL
        if !(fsal isa FSAL)
            ode_wrap!(f_tmp, t + dt/2.0, y2)
            FE[1] += 1
        end
        @.. dy[:,1] = (dt/2.0) * f_tmp
    end
    runge_kutta_step!(method, iteration, y2, t+dt/2.0, dt/2.0, ode_wrap!, dy, y_tmp,
                      f_tmp, FE, error, J, linear_cache, stage_finder)
    @.. y2 = y_tmp
    return nothing
end
