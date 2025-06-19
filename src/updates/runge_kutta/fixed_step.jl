
function evolve_one_time_step!(method::RungeKutta, adaptive::Fixed,
             t::Vector{T}, dt::Vector{T}, ode_wrap_y!::ODEWrapperState,
             update_cache::RKMCache, linear_cache, root_finder::RootFinderMethod,
             stage_finder::ImplicitStageFinder,
             sensitivity::SensitivityMethod, ode_wrap_p!::ODEWrapperParam,
             interpolator::Interpolator) where T <: AbstractFloat

    @unpack iteration, explicit_stage, fesal = method
    @unpack dy, y, y_tmp, f, f_tmp = update_cache

    # evaluate first stage at current time
    if explicit_stage[1]
        @.. dy[:,1] = dt[1] * f
    end

    # for sensitivity
    if explicit_stage[1]
        stage_idx = 1
        @.. y_tmp = y
        explicit_sensitivity_stage!(sensitivity, stage_idx, stage_finder, t[1],
                                    dt[1], update_cache, ode_wrap_y!, ode_wrap_p!,
                                    method)
    end

    runge_kutta_step!(method, iteration, t, dt, ode_wrap_y!, update_cache,
                      linear_cache, root_finder, stage_finder, sensitivity, ode_wrap_p!)

    # evaluate ODE at next time step and store in f_tmp (skip if method is FESAL)
    if (explicit_stage[1] || interpolator isa CubicHermite) && !fesal
        ode_wrap_y!(f_tmp, t[1] + dt[1], y_tmp)
    end

    return nothing
end