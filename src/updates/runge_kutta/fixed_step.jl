
function evolve_one_time_step!(method::RungeKutta, adaptive::Fixed,
             controller::Controller, t::Vector{T}, dt::Vector{T},
             ode_wrap!::ODEWrapperState, update_cache::RKMCache, linear_cache,
             stage_finder::ImplicitStageFinder,
             sensitivity_method::SensitivityMethod,
             ode_wrap_p!::ODEWrapperParam) where T <: AbstractFloat

    @unpack iteration, explicit_stage, fsal = method
    @unpack dy, y, y_tmp, f, f_tmp = update_cache

    # don't want to commit just yet (should check TRBDF2 and Backward Euler)
    if ode_wrap!.FE[1] == 0
        # always evaluate first stage at initial time (should move outside of function)
        ode_wrap!(f, t[1], y)
    else
        # get ODE of current time step (should already be stored in f_tmp)
        @.. f = f_tmp
    end

    @.. dy[:,1] = dt[1] * f

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

    # evaluate ODE at next time step and store in f_tmp (skip if method is FSAL)
    # if (explicit_stage[1] || interpolator isa HermiteInterpolator) && !fsal
    if !fsal
        ode_wrap!(f_tmp, t[1] + dt[1], y_tmp)
    end

    return nothing
end