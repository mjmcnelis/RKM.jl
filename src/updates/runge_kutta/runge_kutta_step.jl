# benchmark.jl note: don't see as much benefit to @muladd as @..
# benchmark.jl note: @.. doesn't help much when switch to MVector
@muladd function runge_kutta_step!(method::RungeKutta, iteration::Explicit,
                                   t::Vector{T}, dt::Vector{T},
                                   config::RKMConfig) where T <: AbstractFloat

    ode_wrap_y! = config.ode_wrap_y!
    update_cache = config.update_cache
    sensitivity = config.sensitivity

    c = method.c
    A_T = method.A_T
    b = method.b
    stages = method.stages

    dy = update_cache.dy
    y = update_cache.y
    y_tmp = update_cache.y_tmp
    f_tmp = update_cache.f_tmp

    for i in 2:stages                                    # evaluate remaining stages
        t_tmp = t[1] + c[i]*dt[1]                        # assumes first stage pre-evaluated
        @.. y_tmp = y
        # TODO: need a better dy cache for performance
        for j in 1:i-1
            if iszero(A_T[j,i])
                continue
            end
            @.. y_tmp = y_tmp + A_T[j,i]*dy[:,j]
        end
        # TODO: skip if intermediate update not needed in next row(s)?
        ode_wrap_y!(f_tmp, t_tmp, y_tmp)
        @.. dy[:,i] = dt[1] * f_tmp

        # for sensitivity
        explicit_sensitivity_stage!(sensitivity, i, t_tmp, dt[1], config, method)
    end
    @.. y_tmp = y                                        # evaluate iteration
    for j in 1:stages
        if iszero(b[j])
            continue
        end
        @.. y_tmp = y_tmp + b[j]*dy[:,j]
    end
    # TODO: just make a function
    if !(sensitivity isa NoSensitivity)
        S = update_cache.S
        S_tmp = update_cache.S_tmp
        dS = update_cache.dS

        @.. S_tmp = S
        for j in 1:stages
            if iszero(b[j])
                continue
            end
            @.. S_tmp = S_tmp + b[j]*dS[:,:,j]
        end
    end
    return nothing
end

@muladd function runge_kutta_step!(method::RungeKutta, iteration::DiagonalImplicit,
                                   t::Vector{T}, dt::Vector{T},
                                   config::RKMConfig) where T <: AbstractFloat

    ode_wrap_y! = config.ode_wrap_y!
    update_cache = config.update_cache
    state_jacobian = config.state_jacobian
    root_finder = config.root_finder
    eigenmax = config.eigenmax
    sensitivity = config.sensitivity

    c = method.c
    A_T = method.A_T
    b = method.b
    stages = method.stages
    explicit_stage = method.explicit_stage

    epsilon = root_finder.epsilon
    p_norm = root_finder.p_norm
    max_iterations = root_finder.max_iterations

    dy = update_cache.dy
    y = update_cache.y
    y_tmp = update_cache.y_tmp
    f_tmp = update_cache.f_tmp
    J = update_cache.J
    res = update_cache.res
    lambda_LR = update_cache.lambda_LR
    x0 = update_cache.x0

    # have while loop stashed but couldn't figure out why it was allocating
    #= while loop: can try something like
        restart_step = true
        while restart_step
            ...
            restart_step = false    # if it finishes
        end
    =#
    # so far some of the changes I put in by hand seem fine

    for i in 1:stages
        # first explicit stage should already be pre-evaluated elsewhere
        if i == 1 && explicit_stage[i]
            # note: redo this if restart step with new dt (attempt > 0)
            # TODO: can probably skip this the first time (attempt = 0)
            #=
            attempt = 0     # TMP
            if attempt > 0
                @.. dy[:,i] = dt[1] * f
                @.. y_tmp = y
                explicit_sensitivity_stage!(sensitivity, i, state_jacobian, t[1], dt[1],
                                            update_cache, ode_wrap_y!, ode_wrap_p!, method)
            end
            =#
            continue
        end

        # set intermediate time in wrapper
        t_tmp = t[1] + c[i]*dt[1]
        set_wrapper!(ode_wrap_y!, t_tmp)

        if explicit_stage[i]
            @.. y_tmp = y
             for j in 1:i-1
                if iszero(A_T[j,i])
                    continue
                end
                @.. y_tmp = y_tmp + A_T[j,i]*dy[:,j]
            end
            ode_wrap_y!(f_tmp, t_tmp, y_tmp)
            @.. dy[:,i] = dt[1] * f_tmp
        else
            # trivial predictor
            @.. dy[:,i] = 0.0

            for n in 1:max_iterations+1
                # evaluate current correction and ODE
                @.. y_tmp = y
                for j in 1:i
                    if iszero(A_T[j,i])
                        continue
                    end
                    @.. y_tmp = y_tmp + A_T[j,i]*dy[:,j]
                end
                ode_wrap_y!(f_tmp, t_tmp, y_tmp)

                # compute residual error of root equation
                # dy - dt.f(t_tmp, y_tmp + A.dy) = 0
                @.. res = dy[:,i] - dt[1]*f_tmp

                # check for convergence (after at least one iteration)
                if n > 1
                    # compute norms and tolerance
                    # note: LinearAlgebra.norm is slow on mac (so use AppleAccelerate)
                    e_norm = norm(res, p_norm)
                    dy_norm = norm(view(dy,:,i), p_norm)
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

                # TODO: if J is sparse, make sure all diagonals are nonzero
                if root_finder isa Newton
                    # evaluate current state Jacobian
                    evaluate_jacobian!(state_jacobian, J, ode_wrap_y!, y_tmp, f_tmp)
                    # note: might want to compute eigenvalues before this

                    # J <- I - dt*A*J
                    root_jacobian!(J, A_T[i,i], dt[1])
                end

                root_iteration!(root_finder, dy, i, res, J)
            end
        end

        # for sensitivity
        implicit_sensitivity_stage!(sensitivity, i, t_tmp, dt[1], config, A_T[i,i], method)
    end
    # evaluate update
    @.. y_tmp = y
    for j in 1:stages
        if iszero(b[j])
            continue
        end
        @.. y_tmp = y_tmp + b[j]*dy[:,j]
    end
    # TODO: just make a function
    if !(sensitivity isa NoSensitivity)
        S = update_cache.S
        S_tmp = update_cache.S_tmp
        dS = update_cache.dS

        @.. S_tmp = S
        for j in 1:stages
            if iszero(b[j])
                continue
            end
            @.. S_tmp = S_tmp + b[j]*dS[:,:,j]
        end
    end

    # estimate max eigenvalue of jacobian
    compute_max_eigenvalue!(eigenmax, lambda_LR, x0, J, state_jacobian,
                            ode_wrap_y!, y_tmp, f_tmp)

    return nothing
end