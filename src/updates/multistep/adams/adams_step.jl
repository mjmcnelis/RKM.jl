
@muladd function adams_step!(method::Adams, ::Explicit,
                     t::T, dt::T, ode_wrap!::ODEWrapperState,
                     update_cache::RKMCache, args...) where T <: AbstractFloat

    @unpack b, stages = method
    @unpack dy_LM, y, y_tmp = update_cache

    # evaluate update
    @.. y_tmp = y
    for j in 1:stages
        dy_stage = view(dy_LM,:,j)
        @.. y_tmp = y_tmp + b[j]*dy_stage
    end
    return nothing
end

@muladd function adams_step!(method::Adams, ::SingleImplicit,
                     t::T, dt::T, ode_wrap!::ODEWrapperState,
                     update_cache::RKMCache, linear_cache,
                     stage_finder::ImplicitStageFinder) where T <: AbstractFloat

    @unpack b, b_pred, stages = method
    @unpack root_method, state_jacobian, epsilon, max_iterations, p_norm = stage_finder
    @unpack dy_LM, y, y_tmp, f_tmp, J, res = update_cache

    # set implicit time in wrapper
    t_tmp = t + dt
    set_wrapper!(ode_wrap!, t_tmp)

    # compute predictor and evaluate ODE (i.e. PE)
    @.. y_tmp = y
    # note: comment for loop if want to compare to BackwardEuler1, ImplicitTrapezoid21
    for j in 1:stages
        dy_stage = view(dy_LM,:,j)
        @.. y_tmp = y_tmp + b_pred[j]*dy_stage
    end
    ode_wrap!(f_tmp, t_tmp, y_tmp)
    @.. dy_LM[:,1] = dt * f_tmp

    for n in 1:max_iterations+1
        # compute current correction and evaluate ODE (i.e CE)
        @.. y_tmp = y
        for j in 1:stages
            dy_stage = view(dy_LM,:,j)
            @.. y_tmp = y_tmp + b[j]*dy_stage
        end
        ode_wrap!(f_tmp, t_tmp, y_tmp)

        # compute residual error of root equation
        # dy - dt.f(t_tmp, y_tmp + b.dy) = 0
        dy_stage = view(dy_LM,:,1)
        @.. res = dy_stage - dt*f_tmp

        # compute norms and tolerance
        e_norm  = norm(res, p_norm)           # compute norms
        dy_norm = norm(dy_stage, p_norm)
        tol = epsilon * dy_norm

        # check for root convergence
        if e_norm < tol
            break
        elseif n == max_iterations + 1
            # println("failed to converge after $(n-1) iteration(s)")
            break
        end

        if root_method isa FixedPoint
            @.. dy_LM[:,1] -= res
        elseif root_method isa Newton
            # evaluate current Jacobian
            evaluate_jacobian!(state_jacobian, J, ode_wrap!, y_tmp, f_tmp)
            # J <- I - b.dt.J
            @.. J *= (-b[1]*dt)
            for k in diagind(J)
                J[k] = J[k] + 1.0
            end

            # pass Jacobian and residual error to linear cache
            linear_cache.A = J
            linear_cache.b = res
            # solve and apply Newton iteration
            solve!(linear_cache)
            @.. dy_LM[:,1] -= linear_cache.u
        end
    end
    # evaluate update
    @.. y_tmp = y
    for j in 1:stages
        dy_stage = view(dy_LM,:,j)
        @.. y_tmp = y_tmp + b[j]*dy_stage
    end
    return nothing
end