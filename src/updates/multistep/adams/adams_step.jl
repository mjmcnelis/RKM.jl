
@muladd function adams_step!(method::Adams, ::Explicit,
                     t::T, dt::T, ode_wrap!::ODEWrapperState,
                     update_cache::RKMCache, args...) where T <: AbstractFloat

    b = method.b
    stages = method.stages

    dy_LM = update_cache.dy_LM
    y = update_cache.y
    y_tmp = update_cache.y_tmp

    # evaluate update
    @.. y_tmp = y
    for j in 1:stages
        @.. y_tmp = y_tmp + b[j]*dy_LM[:,j]
    end
    return nothing
end

@muladd function adams_step!(method::Adams, ::SingleImplicit,
                     t::T, dt::T, ode_wrap!::ODEWrapperState,
                     update_cache::RKMCache, linear_cache,
                     state_jacobian::JacobianMethod, root_finder::RootFinderMethod,
                     eigenmax::EigenMaxMethod) where T <: AbstractFloat

    b = method.b
    b_pred = method.b_pred
    stages = method.stages

    epsilon = root_finder.epsilon
    p_norm = root_finder.p_norm
    max_iterations = root_finder.max_iterations

    dy_LM = update_cache.dy_LM
    y = update_cache.y
    y_tmp = update_cache.y_tmp
    f_tmp = update_cache.f_tmp
    J = update_cache.J
    res = update_cache.res

    # set implicit time in wrapper
    t_tmp = t + dt
    set_wrapper!(ode_wrap!, t_tmp)

    # compute predictor and evaluate ODE (i.e. PE)
    @.. y_tmp = y
    # note: comment for loop if want to compare to BackwardEuler1, ImplicitTrapezoid21
    for j in 1:stages
        @.. y_tmp = y_tmp + b_pred[j]*dy_LM[:,j]
    end
    ode_wrap!(f_tmp, t_tmp, y_tmp)
    @.. dy_LM[:,1] = dt * f_tmp

    for n in 1:max_iterations+1
        # compute current correction and evaluate ODE (i.e CE)
        @.. y_tmp = y
        for j in 1:stages
            @.. y_tmp = y_tmp + b[j]*dy_LM[:,j]
        end
        ode_wrap!(f_tmp, t_tmp, y_tmp)

        # compute residual error of root equation
        # dy - dt.f(t_tmp, y_tmp + b.dy) = 0
        @.. res = dy_LM[:,1] - dt*f_tmp

        # compute norms and tolerance
        e_norm  = norm(res, p_norm)           # compute norms
        dy_norm = norm(view(dy_LM,:,1), p_norm)
        tol = epsilon * dy_norm

        # check for root convergence
        if e_norm < tol
            break
        elseif n == max_iterations + 1
            # println("failed to converge after $(n-1) iteration(s)")
            break
        end

        if root_finder isa FixedPoint
            @.. dy_LM[:,1] -= res
        elseif root_finder isa Newton
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
        @.. y_tmp = y_tmp + b[j]*dy_LM[:,j]
    end
    return nothing
end