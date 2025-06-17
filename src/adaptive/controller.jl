
function initialize_controller!(update_cache::UpdateCache, e_norm::T,
                                tol::T, dt::T) where T <: AbstractFloat

    @unpack e_prev, tol_prev, dt_prev = update_cache

    e_prev .= e_norm
    tol_prev .= tol
    dt_prev .= dt

    return nothing
end

function rescale_time_step(adaptive::ATS, update_cache::UpdateCache, tol::T,
                           e_norm::T) where {ATS <: AdaptiveTimeStep, T <: AbstractFloat}

    @unpack e_prev, tol_prev, dt_prev = update_cache
    @unpack pid = adaptive
    @unpack beta1, beta2, beta3, alpha2, alpha3 = pid

    rescale = (tol/e_norm)^beta1 * (tol_prev[1]/e_prev[1])^beta2 *
              (tol_prev[2]/e_prev[2])^beta3 * (dt_prev[2]/dt_prev[1])^alpha2 *
              (dt_prev[3]/dt_prev[2])^alpha3

    return rescale
end

function set_previous_control_vars!(update_cache::UpdateCache, e_norm::T,
                                    tol::T, dt::T) where T <: AbstractFloat

    @unpack e_prev, tol_prev, dt_prev = update_cache

    e_prev[2]   = e_prev[1]
    e_prev[1]   = e_norm
    tol_prev[2] = tol_prev[1]
    tol_prev[1] = tol
    dt_prev[3]  = dt_prev[2]
    dt_prev[2]  = dt_prev[1]
    dt_prev[1]  = dt

    return nothing
end

function adjust_final_time_steps!(t::Vector{T}, dt::Vector{T},
                                  tf::T) where T <: AbstractFloat
    if dt[2] > tf - t[1]
        dt[2] = tf - t[1]
        dt[1] = dt[2]
    end
    return nothing
end