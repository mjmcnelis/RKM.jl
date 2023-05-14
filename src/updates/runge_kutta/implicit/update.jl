
@muladd function fixed_runge_kutta_step!(method::RungeKutta, ::DiagonalImplicit,
                     y::VectorMVector, t::T, dt::T, dy_dt!::F, dy::MatrixMMatrix,
                     y_tmp::VectorMVector, f_tmp::VectorMVector, J::MatrixMMatrix,
                     linear_cache, finitediff_cache, jacobian_config, 
                     dy_dt_wrap!) where {T <: AbstractFloat, F <: Function}

    @unpack c, A_T, b, stages = method

    # root_solver = "fixed_point"        # will just use fixed point iteration for now
    root_solver = "newton"
    # root_solver = "newton_fast"
    eps_root = 1e-8
    max_iterations = 2                   # 2 for benchmarking ImplicitEuler(), fixed dt

    for i = 1:stages
        # evaluate jacobian
        #=
        if root_solver == "newton_fast" 
            # TODO: how to reduce allocations here, take it apart or make a wrapper?
            # so in order for jacobian config to work, dy_dt_wrap! argument
            # has to be the same object stored in jacobian_config
            jacobian!(J, dy_dt_wrap!, f_tmp, y, jacobian_config)
            # @show J
            # finite_difference_jacobian!(J, dy_dt_wrap!, y, finitediff_cache) 
            # @show J; q()

            J .*= (-A_T[i,i]*dt)
            for i in diagind(J)
                J[i] += 1.0
            end
            # pass jacobian to linear cache
            linear_cache = set_A(linear_cache, J)
        end
        =#
        t_tmp = t + c[i]*dt
        # sum over previously known stages
        @.. y_tmp = y 
        for j = 1:i-1
            dy_stage = view(dy,:,j)
            @.. y_tmp = y_tmp + A_T[j,i]*dy_stage
        end
        # zero current stage before iterating
        @.. dy[:,i] = 0.0

        # TEMP iterate w/o any breaks for now
        for n = 1:max_iterations
            dy_stage = view(dy,:,i)
            @.. y_tmp = y_tmp + A_T[i,i]*dy_stage

            # do this before or after y_tmp evaluation? 
            if root_solver == "newton"
                jacobian!(J, dy_dt_wrap!, f_tmp, y_tmp, jacobian_config)
                J .*= (-A_T[i,i]*dt)
                for i in diagind(J)
                    J[i] += 1.0
                end
                linear_cache = set_A(linear_cache, J)
            end
        
            # evaluate ODE at previous iteration (could maybe just iterate y_tmp itself)
            dy_dt!(f_tmp, t_tmp, y_tmp)

            # TEMP undo addition (for minimizing allocations)
            @.. y_tmp = y_tmp - A_T[i,i]*dy_stage

            if root_solver == "fixed_point"
                @.. dy[:,i] = dt * f_tmp
            elseif root_solver == "newton"
                # TODO: try to solve for dy directly instead of d(dy)
                for k in eachindex(f_tmp)
                    # TODO: shouldn't this include the stage coefficient?
                    #       maybe not, double check stage math
                    f_tmp[k] = dy[k,i] - dt*f_tmp[k]
                end
                # from python
                # g = z - dt*y_prime(t + dt*c[i], y + dy + z*Aii)

                linear_cache = set_b(linear_cache, f_tmp)
                # note: may not need this if use regular newton method 
                linear_cache = solve_linear_tmp(linear_cache)
                @.. dy[:,i] -= linear_cache.u
                # sol = solve(linear_cache)
                # @.. dy[:,i] -= sol.u
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
            J::MatrixMMatrix, linear_cache, finitediff_cache,
            jacobian_config, dy_dt_wrap!, args...) where {T <: AbstractFloat, F}

    fixed_runge_kutta_step!(method, iteration, y, t[1], dt[1], dy_dt!, dy, y_tmp, f_tmp,
                            J, linear_cache, finitediff_cache, jacobian_config,
                            dy_dt_wrap!)
    y .= y_tmp
    nothing
end
function evolve_one_time_step!(method::RungeKutta, iteration::DiagonalImplicit,
             adaptive::Doubling, controller::Controller, FE::MVector{1,Int64},  
             y::VectorMVector, t::VectorMVector{1,T}, dt::VectorMVector{2,T},
             dy_dt!::F, dy::MatrixMMatrix, y_tmp::VectorMVector, 
             f_tmp::VectorMVector, f::VectorMVector, y1::VectorMVector, 
             y2::VectorMVector, error::VectorMVector, 
             J::MatrixMMatrix, linear_cache, finitediff_cache, 
             jacobian_config, dy_dt_wrap!, args...) where {T <: AbstractFloat, F}

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

        doubling_runge_kutta_step!(method, iteration, y, t[1], dt[1],
                                dy_dt!, dy, y_tmp, f_tmp, f, y1, y2, J, 
                                linear_cache, finitediff_cache, jacobian_config, 
                                dy_dt_wrap!)

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

function doubling_runge_kutta_step!(method, iteration::DiagonalImplicit, y, t, dt,
                                    dy_dt!, dy, y_tmp, f_tmp, f, y1, y2, J, 
                                    linear_cache, finitediff_cache, jacobian_config, 
                                    dy_dt_wrap!)
    # iterate full time step 
    fixed_runge_kutta_step!(method, iteration, y, t, dt, dy_dt!, dy, y_tmp, f_tmp, 
                            J, linear_cache, finitediff_cache, jacobian_config, 
                            dy_dt_wrap!)
    @.. y1 = y_tmp                                      # y1(t+dt)   
    # iterate two half time steps
    fixed_runge_kutta_step!(method, iteration, y, t, dt/2.0, dy_dt!, dy, y_tmp, f_tmp,
                            J, linear_cache, finitediff_cache, jacobian_config, 
                            dy_dt_wrap!)
    @.. y2 = y_tmp                                      # y2(t+dt/2)
    fixed_runge_kutta_step!(method, iteration, y2, t + dt/2.0, dt/2.0, dy_dt!, dy, y_tmp, 
                            f_tmp, J, linear_cache, finitediff_cache, jacobian_config, 
                            dy_dt_wrap!)
    @.. y2 = y_tmp                                      # y2(t+dt)
    return nothing
end