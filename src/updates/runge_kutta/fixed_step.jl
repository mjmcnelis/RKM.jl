
function evolve_one_time_step!(method::RungeKutta, adaptive::Fixed,
             controller::Controller, t::Vector{T}, dt::Vector{T},
             ode_wrap!::ODEWrapperState, update_cache::RKMCache, linear_cache,
             stage_finder::ImplicitStageFinder,
             sensitivity_method::SensitivityMethod,
             ode_wrap_p!::ODEWrapperParam) where T <: AbstractFloat

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
        explicit_sensitivity_stage!(sensitivity_method, stage_idx, stage_finder,
                                    t[1], dt[1], update_cache, ode_wrap!,
                                    ode_wrap_p!)
    end

    runge_kutta_step!(method, iteration, t[1], dt[1], ode_wrap!,
                      update_cache, linear_cache, stage_finder,
                      sensitivity_method, ode_wrap_p!)

    # evaluate ODE at next time step and store in f_tmp (skip if method is FESAL)
    # note: get excess allocations if try to pass interpolator
    # if (explicit_stage[1] || interpolator isa HermiteInterpolator) && !fesal
    if !fesal
        ode_wrap!(f_tmp, t[1] + dt[1], y_tmp)
    end

    return nothing
end