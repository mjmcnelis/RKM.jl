
function evolve_one_time_step!(method::MagnusAugmented, adaptive::Fixed,
                               t::Vector{T}, dt::Vector{T},
                               config::RKMConfig) where T <: AbstractFloat

    ode_wrap_y! = config.ode_wrap_y!
    update_cache = config.update_cache
    interpolator = config.interpolator

    magnus_step!(method, t, dt, config)

    # TODO: check if Magnus really needs this
    if interpolator isa CubicHermite
        y_tmp = update_cache.y_tmp
        f_tmp = update_cache.f_tmp
        ode_wrap_y!(f_tmp, t[1] + dt[1], y_tmp)
    end

    return nothing
end