
@muladd function adams_step!(method::Adams, iteration::Explicit,
                             t::Vector{T}, dt::Vector{T},
                             config::RKMConfig) where T <: AbstractFloat

    update_cache = config.update_cache

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

@muladd function adams_step!(method::Adams, iteration::SingleImplicit,
                             t::Vector{T}, dt::Vector{T},
                             config::RKMConfig) where T <: AbstractFloat

    ode_wrap! = config.ode_wrap_y!
    update_cache = config.update_cache
    state_jacobian = config.state_jacobian
    root_finder = config.root_finder

    b = method.b
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

    # set intermediate time in wrapper
    t_tmp = t[1] + dt[1]
    set_wrapper!(ode_wrap!, t_tmp)

    # trivial predictor
    @.. dy_LM[:,1] = 0.0

    for n in 1:max_iterations+1
        # compute current correction and evaluate ODE (i.e CE)
        @.. y_tmp = y
        for j in 1:stages
            @.. y_tmp = y_tmp + b[j]*dy_LM[:,j]
        end
        ode_wrap!(f_tmp, t_tmp, y_tmp)

        # compute residual error of root equation
        # dy - dt.f(t_tmp, y_tmp + b.dy) = 0
        @.. res = dy_LM[:,1] - dt[1]*f_tmp

        # check for convergence (after at least one iteration)
        if n > 1
            # compute norms and tolerance
            # note: LinearAlgebra.norm is slow on mac (so use AppleAccelerate)
            e_norm  = norm(res, p_norm)
            dy_norm = norm(view(dy_LM,:,1), p_norm)
            tol = epsilon * dy_norm

            if e_norm <= tol
                break
            end
        end

        if n == max_iterations + 1
            # println("failed to converge after $(n-1) iteration(s)")
            break
        end

        if root_finder isa Newton
            # evaluate current Jacobian
            evaluate_jacobian!(state_jacobian, J, ode_wrap!, y_tmp, f_tmp)
            # J <- I - b.dt.J
            root_jacobian!(J, b[1], dt[1])
        end

        root_iteration!(root_finder, dy_LM, 1, res, J)
    end

    # evaluate update
    @.. y_tmp = y
    for j in 1:stages
        @.. y_tmp = y_tmp + b[j]*dy_LM[:,j]
    end

    return nothing
end