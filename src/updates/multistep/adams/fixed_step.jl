
function evolve_one_time_step!(method::Adams, adaptive::Fixed,
             controller::Controller, FE::MVector{1,Int64},
             t::VectorMVector{1,T}, dt::VectorMVector{2,T},
             ode_wrap!::ODEWrapper, update_cache::RKMCache, linear_cache,
             stage_finder::ImplicitStageFinder) where T <: AbstractFloat

    @unpack stages, iteration = method
    @unpack dy, y, y_tmp, f = update_cache

    # evaluate ODE at current state/time (i.e E)
    ode_wrap!(f, t[1], y)
    @.. dy[:,1] = dt[1] * f

    # initialize previous stages (trivial)
    if FE[1] == 0
        for j in 2:stages
            @.. dy[:,j] = dt[1] * f
        end
    end
    FE[1] += 1

    adams_step!(method, iteration, t[1], dt[1], ode_wrap!, FE,
                update_cache, linear_cache, stage_finder)

    # shift stages (except first one)
    for j in 1:stages-1
        @.. dy[:,stages+1-j] = dy[:,stages-j]
    end
    @.. y = y_tmp
    return nothing
end