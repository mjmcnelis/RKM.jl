
@muladd function fixed_runge_kutta_step!(method::RungeKutta, ::DiagonalImplicit,
                     y::VectorMVector, t::T, dt::T, dy_dt!::F, dy::MatrixMMatrix,
                     y_tmp::VectorMVector, f_tmp::VectorMVector, J::MatrixMMatrix,
                     linear_cache, dy_dt_wrap!,
                     stage_finder::ImplicitStageFinder) where {T <: AbstractFloat, 
                                                               F <: Function}

    @unpack c, A_T, b, stages, explicit_stage = method
    @unpack root_method, jacobian_method, epsilon, max_iterations = stage_finder

    for i = 1:stages        
        t_tmp = t + c[i]*dt
        @.. y_tmp = y                                    # sum over known stages
        for j = 1:i-1
            dy_stage = view(dy,:,j)
            @.. y_tmp = y_tmp + A_T[j,i]*dy_stage
        end
        
        if explicit_stage[i]
            dy_dt!(f_tmp, t_tmp, y_tmp)
            @.. dy[:,i] = dt * f_tmp
        else
            # TODO: look into predictors
            dy_dt!(f_tmp, t_tmp, y_tmp)                  # guess stage before iterating
            @.. dy[:,i] = dt * f_tmp

            for n = 1:max_iterations
                dy_stage = view(dy,:,i)
                @.. y_tmp = y_tmp + A_T[i,i]*dy_stage

                if root_method isa Newton                # evaluate current Jacobian
                    evaluate_system_jacobian!(jacobian_method, J, dy_dt_wrap!, 
                                              y_tmp, f_tmp)
                    J .*= (-A_T[i,i]*dt)                 # J <- I - A.dt.J
                    for i in diagind(J)
                        J[i] += 1.0
                    end
                    linear_cache = set_A(linear_cache, J)# pass Jacobian to linear cache
                end

                dy_dt!(f_tmp, t_tmp, y_tmp)              # evaluate current slope

                @.. y_tmp = y_tmp - A_T[i,i]*dy_stage    # undo addition to y_tmp

                if root_method isa FixedPoint
                    @.. dy[:,i] = dt * f_tmp
                elseif root_method isa Newton
                    # TODO: initialize predictor for dy other than zero,
                    #       and sort out @.. equivalent
                    for k in eachindex(f_tmp)
                        f_tmp[k] = dy[k,i] - dt*f_tmp[k]
                    end

                    # compute residual error of root equation:
                    # dy - dt.f(t_tmp, y_tmp + A.dy) = 0
                    res = norm(f_tmp)

                    dy_norm = norm(view(dy,:,i))         # compute error tolerance 
                    tol = epsilon * dy_norm
               
                    if n > 1 && res < tol                # check if Newton method covnerges
                        # print(n)
                        break
                    end
                    # if n == max_iterations
                    #     @warn "exceeded max Newton iterations at t = $t"
                    #     @show res tol
                    #     println("")
                    # end

                    linear_cache = set_b(linear_cache, f_tmp)
                    # note: may not need this if use regular newton method 
                    linear_cache = solve_linear_tmp(linear_cache)
                    @.. dy[:,i] -= linear_cache.u
                    # sol = solve(linear_cache)
                    # @.. dy[:,i] -= sol.u
                end
            end
        end
    end
    @.. y_tmp = y                                        # evaluate iteration
    for j = 1:stages
        dy_stage = view(dy,:,j)
        @.. y_tmp = y_tmp + b[j]*dy_stage
    end
    nothing
end

# TODO: use alternative indents
function evolve_one_time_step!(method::RungeKutta, iteration::DiagonalImplicit,
            adaptive::Fixed, ::Controller, FE::MVector{1,Int64},
            y::VectorMVector, t::VectorMVector{1,T}, dt::VectorMVector{2,T},
            dy_dt!::F, dy::MatrixMMatrix, y_tmp::VectorMVector, 
            f_tmp::VectorMVector, f::VectorMVector, y1, y2, error, 
            J::MatrixMMatrix, linear_cache, dy_dt_wrap!,
            stage_finder::ImplicitStageFinder) where {T <: AbstractFloat, F}

    fixed_runge_kutta_step!(method, iteration, y, t[1], dt[1], dy_dt!, dy, y_tmp, 
                            f_tmp, J, linear_cache, dy_dt_wrap!, stage_finder)
    y .= y_tmp
    nothing
end
function evolve_one_time_step!(method::RungeKutta, iteration::DiagonalImplicit,
             adaptive::Doubling, controller::Controller, FE::MVector{1,Int64},  
             y::VectorMVector, t::VectorMVector{1,T}, dt::VectorMVector{2,T},
             dy_dt!::F, dy::MatrixMMatrix, y_tmp::VectorMVector, 
             f_tmp::VectorMVector, f::VectorMVector, y1::VectorMVector, 
             y2::VectorMVector, error::VectorMVector, 
             J::MatrixMMatrix, linear_cache, dy_dt_wrap!, 
             stage_finder::ImplicitStageFinder) where {T <: AbstractFloat, F}

    @unpack epsilon, low, high, safety, p_norm, dt_min, dt_max, max_attempts = adaptive

    order = method.order[1]                             # order of scheme

    # note: comment if want to compare to OrdinaryDiffEq ImplicitEuler()
    high ^= order / (1.0 + order)                       # rescale high based on order
    epsilon ^= order / (1.0 + order)

    dt[1] = dt[2]                                       # initialize time step
    rescale = 1.0

    attempts = 1
    while true                                          # start step doubling routine
        dt[1] = min(dt_max, max(dt_min, dt[1]*rescale)) # increase dt for next attempt

        doubling_runge_kutta_step!(method, iteration, y, t[1], dt[1], dy_dt!, dy, y_tmp,
                                   f_tmp, f, y1, y2, J, linear_cache, dy_dt_wrap!, 
                                   stage_finder)

        @.. error = (y2 - y1) / (2.0^order - 1.0)       # estimate local truncation error
        @.. y2 = y2 + error                             # Richardson extrapolation

        e_norm = norm(error, p_norm)                    # compute norms
        y_norm = norm(y2, p_norm)
        Δy = y1
        @.. Δy = y2 - y
        Δy_norm = norm(Δy, p_norm)

        tol = epsilon * max(y_norm, Δy_norm)            # compute tolerance

        if FE[1] == 0                                   # initialize controller
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
    @.. y = y2                                          # get iteration
    # add_function_evaluations!(FE, iteration, adaptive, method, attempts)
    FE[1] += 10
    return nothing
end

# TODO: combine explicit and implicit routines
function evolve_one_time_step!(method::RungeKutta, iteration::DiagonalImplicit,
             adaptive::Embedded, controller::Controller, FE::MVector{1,Int64},  
             y::VectorMVector, t::VectorMVector{1,T}, dt::VectorMVector{2,T},
             dy_dt!::F, dy::MatrixMMatrix, y_tmp::VectorMVector, 
             f_tmp::VectorMVector, f::VectorMVector, y1::VectorMVector, 
             y2::VectorMVector, error::VectorMVector, 
             J::MatrixMMatrix, linear_cache, dy_dt_wrap!, 
             stage_finder::ImplicitStageFinder) where {T <: AbstractFloat, F}

    @unpack epsilon, low, high, safety, p_norm, dt_min, dt_max, max_attempts = adaptive

    y_norm = norm(y, p_norm)                            # compute norm of current state

    order_max = maximum(method.order)                   # max/min orders in embedded scheme
    order_min = minimum(method.order)

    # note: comment if benchmark against OrdinaryDiffEq (add boolean?)
    high    ^= (order_min / order_max)                  # rescale high, epsilon parameters
    epsilon ^= (order_min / order_max)

    dt[1] = dt[2]                                       # initialize time step
    rescale = 1.0

    # # TODO: dispatch fsal
    # if !(method.fsal isa FSAL && FE[1] > 0)
    #     dy_dt!(f, t[1], y)                              # evaluate first stage at (t,y)
    # endz

    attempts = 1
    while true                                          # start step doubling routine
        dt[1] = min(dt_max, max(dt_min, dt[1]*rescale)) # increase dt for next attempt
                       
        fixed_runge_kutta_step!(method, iteration, y, t[1], dt[1], dy_dt!, dy, y_tmp, 
                                f_tmp, J, linear_cache, dy_dt_wrap!, stage_finder)
        @.. y1 = y_tmp                                  # primary iteration

        embedded_runge_kutta_step!(method, y, dy, y_tmp)
        @.. y2 = y_tmp                                  # embedded iteration

        @.. error = y2 - y1                             # local error of embedded pair

        e_norm = norm(error, p_norm)                    # compute norms
        y1_norm = norm(y1, p_norm)
        # TODO: need to use Δy = y2 since it's secondary method
        #       but should make labeling consistent w/ doubling
        Δy = y2
        @.. Δy = y1 - y
        Δy_norm = norm(Δy, p_norm)

        # compute tolerance
        tol = epsilon * max(y_norm, y1_norm, Δy_norm)

        if FE[1] == 0                                   # initialize controller
            initialize_controller!(controller, e_norm, dt[1])
        end

        if e_norm == 0.0                                # compute scaling factor for dt
            rescale = high
        else
            rescale = rescale_time_step(controller, tol, e_norm, order_min)
            rescale = min(high, max(low, safety*rescale))
        end

        dt[2] = min(dt_max, max(dt_min, dt[1]*rescale)) # projected dt for next iteration

        if e_norm <= tol                                # compare error to tolerance
            set_previous_control_vars!(controller, e_norm, dt[1])
            break 
        end
        attempts <= max_attempts || (@warn "embedded exceeded $max_attempts attempts"; break)
        attempts += 1
    end
    # # TODO: dispatch fsal
    # if method.fsal isa FSAL
    #     @.. f = f_tmp 
    # end
    
    @.. y = y1                                          # get iteration
    # add_function_evaluations!(FE, iteration, adaptive, method, attempts)
    FE[1] += 10
    return nothing
end

function doubling_runge_kutta_step!(method, iteration, y, t, dt, dy_dt!, dy, y_tmp, f_tmp, 
                                    f, y1, y2, J, linear_cache, dy_dt_wrap!, stage_finder)
    # iterate full time step 
    fixed_runge_kutta_step!(method, iteration, y, t, dt, dy_dt!, dy, y_tmp, 
                            f_tmp, J, linear_cache, dy_dt_wrap!, stage_finder)
    @.. y1 = y_tmp                                      # y1(t+dt)   
    # iterate two half time steps
    fixed_runge_kutta_step!(method, iteration, y, t, dt/2.0, dy_dt!, dy, y_tmp, 
                            f_tmp, J, linear_cache, dy_dt_wrap!, stage_finder)
    @.. y2 = y_tmp                                      # y2(t+dt/2)
    fixed_runge_kutta_step!(method, iteration, y2, t + dt/2.0, dt/2.0, dy_dt!, dy, y_tmp, 
                            f_tmp, J, linear_cache, dy_dt_wrap!, stage_finder)
    @.. y2 = y_tmp                                      # y2(t+dt)
    return nothing
end