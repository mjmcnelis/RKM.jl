
function magnus_step!(method::MagnusAugmented,
                      t::Vector{T}, dt::Vector{T},
                      config::RKMConfig) where T <: AbstractFloat

    ode_wrap_y! = config.ode_wrap_y!
    update_cache = config.update_cache
    state_jacobian = config.state_jacobian
    root_finder = config.root_finder

    max_iterations = root_finder.max_iterations
    epsilon = root_finder.epsilon
    convergence_failures = root_finder.convergence_failures
    p_norm = root_finder.p_norm

    y = update_cache.y
    y_tmp = update_cache.y_tmp
    f = update_cache.f
    f_tmp = update_cache.f_tmp
    J = update_cache.J
    J_tmp = update_cache.J_tmp
    res = update_cache.res

    ode_wrap_y!(f, t[1], y)
    set_wrapper!(ode_wrap_y!, t[1])
    evaluate_jacobian!(state_jacobian, J, ode_wrap_y!, y, f)

    # TODO: allocate y_dual (maybe change y to dual)
    # y_dual = Dual.(y, f)
    # t_dual = Dual(t[1], 1.0)
    # f_dual = Dual.(f, 1.0)
    # ode_wrap_y!(f_dual, t_dual, y_dual)
    # dfdt = partials.(f_dual, 1)

    @.. y_tmp = y + dt[1]*f #+ 0.5*dt[1]^2*dfdt

    # reuse f?
    b = f - J*y

    set_wrapper!(ode_wrap_y!, t[1] + dt[1])

    for n in 1:max_iterations+1
        @.. res = y_tmp

        ode_wrap_y!(f_tmp, t[1] + dt[1], y_tmp)
        evaluate_jacobian!(state_jacobian, J_tmp, ode_wrap_y!, y_tmp, f_tmp)

        b_tmp = f_tmp - J_tmp*y_tmp

        b_int = 0.5*dt[1]*(b + b_tmp)
        gamma = b_int - dt[1]^2/12.0*(J*b_tmp - J_tmp*b)

        if n == max_iterations + 1
            println("failed to converge after $(n-1) iterations")
            convergence_failures[1] += 1
            break
        end

        # trapezoid is easiest thing to do now
        # z0 = J_tmp
        # @.. z0 += J
        # @.. z0 *= -0.5*dt[1]
        z0 = -0.5*dt[1]*(J + J_tmp) + dt[1]^2/12.0*(J*J_tmp - J_tmp*J)

        y_mag = exp(-z0)*y + (I - exp(-z0))*(z0\gamma)
        @.. y_tmp = y_mag

        # # TODO: is there a better stopping criteria?
        @.. res -= y_tmp
        e_norm = norm(res, p_norm)
        y_norm = norm(y_tmp, p_norm)
        tol = epsilon * (1.0 + y_norm)

        if e_norm <= tol
            root_finder.iteration_bins[n] += 1
            break
        end
    end

    return nothing
end
