
function evolve_one_time_step!(method::RungeKutta,
             adaptive::Fixed, controller::Controller, FE::MVector{1,Int64},
             y::VectorMVector, t::VectorMVector{1,T}, dt::VectorMVector{2,T},
             ode_wrap!::ODEWrapper, update_cache::RKMCache, linear_cache,
             stage_finder::ImplicitStageFinder) where T <: AbstractFloat

    @unpack iteration, explicit_stage, fsal = method
    @unpack dy, y_tmp, f, f_tmp = update_cache

    # TODO: wrap into a function
    # evaluate first stage at (t,y)
    if explicit_stage[1]
        # skip function evaluation if method is FSAL
        if FE[1] > 0 && fsal isa FSAL
            f .= f_tmp
        else
            ode_wrap!(f, t[1], y)
            FE[1] += 1
        end
        @.. dy[:,1] = dt[1] * f
    end

    runge_kutta_step!(method, iteration, y, t[1], dt[1], ode_wrap!, FE,
                      update_cache, linear_cache, stage_finder)
    @.. y = y_tmp
    return nothing
end