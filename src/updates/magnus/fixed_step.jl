
function evolve_one_time_step!(method::M, adaptive::Fixed,
                               t::Vector{T}, dt::Vector{T},
                               config::RKMConfig) where {M <: Magnus, T <: AbstractFloat}

    ode_wrap_y! = config.ode_wrap_y!
    update_cache = config.update_cache
    state_jacobian = config.state_jacobian
    interpolator = config.interpolator

    y = update_cache.y
    f = update_cache.f
    J = update_cache.J

    # initialize ODE and Jacobian
    ode_wrap_y!(f, t[1], y)
    set_wrapper!(ode_wrap_y!, t[1])
    evaluate_jacobian!(state_jacobian, J, ode_wrap_y!, y, f)

    magnus_step!(method, t, dt, config)

    # TODO: check if Magnus really needs this
    if interpolator isa CubicHermite
        y_tmp = update_cache.y_tmp
        f_tmp = update_cache.f_tmp
        ode_wrap_y!(f_tmp, t[1] + dt[1], y_tmp)
    end

    return nothing
end