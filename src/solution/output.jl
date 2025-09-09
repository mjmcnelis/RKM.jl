
function _output_solution!(sol::Solution{T}, update_cache::UpdateCache{T}, t::Vector{T},
                           options::SolverOptions{T}) where T <: AbstractFloat

    save_time_derivative = options.save_time_derivative
    interpolator = options.interpolator
    sensitivity = options.sensitivity

    push!(sol.t, t[1])
    append!(sol.y, update_cache.y_tmp)
    if save_time_derivative || interpolator isa CubicHermite
        append!(sol.f, update_cache.f_tmp)
    end
    if interpolator isa ContinuousFormula
        append!(sol.dy, update_cache.dy)
    end
    if !(sensitivity isa NoSensitivity)
        append!(sol.S, update_cache.S_tmp)
    end
    return nothing
end

function output_solution!(sol::Solution{T}, save_time::Vector{Float64},
                          update_cache::UpdateCache{T}, t::Vector{T},
                          options::SolverOptions{T},
                          timer::TimeLimit) where T <: AbstractFloat

    time_subroutine = options.time_subroutine
    total_steps = timer.total_steps

    if time_subroutine && total_steps[1] % SAMPLE_INTERVAL == 0
        stats = @timed _output_solution!(sol, update_cache, t, options)
        save_time[1] += SAMPLE_INTERVAL*stats.time
    else
        _output_solution!(sol, update_cache, t, options)
    end
    return nothing
end