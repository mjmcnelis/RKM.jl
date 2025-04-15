# benchmark.jl note: don't see as much benefit to @muladd as @..
# benchmark.jl note: @.. doesn't help much when switch to MVector
@muladd function runge_kutta_step!(method::RungeKutta, ::Explicit,
                     t::T, dt::T, ode_wrap_y!::ODEWrapperState,
                     update_cache::RKMCache, linear_cache,
                     stage_finder::ImplicitStageFinder,
                     sensitivity::SensitivityMethod,
                     ode_wrap_p!::ODEWrapperParam) where T <: AbstractFloat
    @unpack c, A_T, b, stages = method
    @unpack dy, y, y_tmp, f_tmp, S, S_tmp, dS = update_cache

    for i in 2:stages                                    # evaluate remaining stages
        t_tmp = t + c[i]*dt                             # assumes first stage pre-evaluated
        @.. y_tmp = y
        # TODO: need a better dy cache for performance
        for j in 1:i-1
            A_T[j,i] == 0.0 ? continue : nothing
            dy_stage = view(dy,:,j)
            @.. y_tmp = y_tmp + A_T[j,i]*dy_stage
        end
        # TODO: skip if intermediate update not needed in next row(s)?
        ode_wrap_y!(f_tmp, t_tmp, y_tmp)
        @.. dy[:,i] = dt * f_tmp

        # for sensitivity
        explicit_sensitivity_stage!(sensitivity, i, stage_finder, t_tmp, dt,
                                    update_cache, ode_wrap_y!, ode_wrap_p!,
                                    method)
    end
    @.. y_tmp = y                                        # evaluate iteration
    for j in 1:stages
        b[j] == 0.0 ? continue : nothing
        dy_stage = view(dy,:,j)
        @.. y_tmp = y_tmp + b[j]*dy_stage
    end
    # TODO: just make a function
    if !(sensitivity isa NoSensitivity)
        @.. S_tmp = S
        for j in 1:stages
            b[j] == 0.0 ? continue : nothing
            dS_stage = view(dS,:,:,j)
            @.. S_tmp = S_tmp + b[j]*dS_stage
        end
    end
    return nothing
end

@muladd function runge_kutta_step!(method::RungeKutta, ::DiagonalImplicit,
                     t::T, dt::T, ode_wrap_y!::ODEWrapperState,
                     update_cache::RKMCache, linear_cache,
                     stage_finder::ImplicitStageFinder,
                     sensitivity::SensitivityMethod,
                     ode_wrap_p!::ODEWrapperParam) where T <: AbstractFloat

    @unpack c, A_T, b, stages, explicit_stage, fesal = method
    @unpack root_method, state_jacobian, epsilon, max_iterations, p_norm = stage_finder
    @unpack dy, y, y_tmp, f_tmp, J, res, S, S_tmp, dS = update_cache

    for i in 1:stages
        # first explicit stage should already be pre-evaluated elsewhere
        if i == 1 && explicit_stage[i]
            continue
        end

        # set intermediate time in wrapper
        t_tmp = t + c[i]*dt
        set_wrapper!(ode_wrap_y!, t_tmp)

        # sum over known stages
        @.. y_tmp = y
        for j in 1:i-1
            dy_stage = view(dy,:,j)
            @.. y_tmp = y_tmp + A_T[j,i]*dy_stage
        end

        # guess stage before iterating
        # TODO: look into predictors
        ode_wrap_y!(f_tmp, t_tmp, y_tmp)
        @.. dy[:,i] = dt * f_tmp

        # get diagonal coefficient
        A = A_T[i,i]

        if !explicit_stage[i]
            for n in 1:max_iterations+1
                # evaluate current correction and ODE
                @.. y_tmp = y
                for j in 1:i
                    dy_stage = view(dy,:,j)
                    @.. y_tmp = y_tmp + A_T[j,i]*dy_stage
                end
                ode_wrap_y!(f_tmp, t_tmp, y_tmp)

                # compute residual error of root equation
                # dy - dt.f(t_tmp, y_tmp + A.dy) = 0
                dy_stage = view(dy,:,i)
                @.. res = dy_stage - dt*f_tmp

                # check for convergence (after at least one iteration)
                if n > 1
                    # compute norms and tolerance
                    # note: LinearAlgebra.norm is slow on mac (so use AppleAccelerate)
                    e_norm = norm(res, p_norm)
                    dy_norm = norm(dy_stage, p_norm)
                    tol = epsilon * dy_norm

                    if e_norm <= tol
                        break
                    end
                end

                if n == max_iterations + 1
                    # println("failed to converge after $(n-1) iterations")
                    # note: allow f_tmp to store dy_dt! of last iteration for FESAL methods
                    # TODO: count convergence failures in stats
                    break
                end

                if root_method isa FixedPoint
                    @.. dy[:,i] -= res
                elseif root_method isa Newton
                    # evaluate current Jacobian
                    evaluate_jacobian!(state_jacobian, J, ode_wrap_y!, y_tmp, f_tmp)

                    # TODO: make rescale methods
                    if J isa SparseMatrixCSC                # J <- I - A.dt.J
                        @.. J.nzval = J.nzval * (-A*dt)
                    else
                        # note: allocates if J is sparse
                        @.. J = J * (-A*dt)
                    end

                    for k in diagind(J)
                        J[k] = J[k] + 1.0
                    end

                    # pass Jacobian and residual error to linear cache
                    linear_cache.A = J
                    linear_cache.b = res
                    solve!(linear_cache)
                    @.. dy[:,i] -= linear_cache.u
                end
            end
        end

        # for sensitivity
        implicit_sensitivity_stage!(sensitivity, i, stage_finder, t_tmp, dt,
                                    update_cache, ode_wrap_y!, ode_wrap_p!, A,
                                    method)
    end
    # evaluate update
    @.. y_tmp = y
    for j in 1:stages
        dy_stage = view(dy,:,j)
        @.. y_tmp = y_tmp + b[j]*dy_stage
    end
    # TODO: just make a function
    if !(sensitivity isa NoSensitivity)
        @.. S_tmp = S
        for j in 1:stages
            b[j] == 0.0 ? continue : nothing
            dS_stage = view(dS,:,:,j)
            @.. S_tmp = S_tmp + b[j]*dS_stage
        end
    end
    return nothing
end