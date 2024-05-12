
function evolve_one_time_step!(method::RungeKutta, adaptive::Fixed,
             controller::Controller, FE::MVector{1,Int64},
             t::Vector{T}, dt::Vector{T},
             ode_wrap!::ODEWrapper, update_cache::RKMCache, linear_cache,
             stage_finder::ImplicitStageFinder) where T <: AbstractFloat

    @unpack iteration, explicit_stage, fsal = method
    @unpack dy, y, y_tmp, f, f_tmp = update_cache

    # don't want to commit just yet (should check TRBDF2 and Backward Euler)

    if FE[1] == 0
        # always evaluate first stage at initial time (should move outside of function)
        ode_wrap!(f, t[1], y)
        FE[1] += 1
    else
        # get ODE of current time step (should already be stored in f_tmp)
        @.. f = f_tmp
    end

    @.. dy[:,1] = dt[1] * f

    runge_kutta_step!(method, iteration, t[1], dt[1], ode_wrap!, FE,
                      update_cache, linear_cache, stage_finder)

    # evaluate ODE at next time step and store in f_tmp (skip if method is FSAL)
    # if (explicit_stage[1] || interpolator isa HermiteInterpolator) && !fsal
    if !fsal
        ode_wrap!(f_tmp, t[1] + dt[1], y_tmp)
        FE[1] += 1
    end

    return nothing
end