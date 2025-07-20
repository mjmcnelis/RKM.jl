
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

    benchmark_subroutines = options.benchmark_subroutines
    total_steps = timer.total_steps

    if benchmark_subroutines && total_steps[1] % 10 == 0
        save_stat = @timed begin
            _output_solution!(sol, update_cache, t, options)
        end
        save_time[1] += 10*save_stat.time
    else
        _output_solution!(sol, update_cache, t, options)
    end
    return nothing
end