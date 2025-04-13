
function _output_solution!(sol::Solution{T}, update_cache::UpdateCache{T}, t::Vector{T},
                           options::SolverOptions{T}) where T <: AbstractFloat

    @unpack save_time_derivative, interpolator, sensitivity = options
    @unpack y_tmp, f_tmp, dy, S_tmp = update_cache

    append!(sol.t, t[1])
    append!(sol.y, y_tmp)
    if save_time_derivative || interpolator isa CubicHermite
        append!(sol.f, f_tmp)
    end
    if interpolator isa ContinuousFormula
        append!(sol.dy, dy)
    end
    if !(sensitivity isa NoSensitivity)
        append!(sol.S, S_tmp)
    end
    return nothing
end

function output_solution!(sol::Solution{T}, save_time::Vector{Float64},
                          update_cache::UpdateCache{T}, t::Vector{T},
                          options::SolverOptions{T},
                          timer::TimeLimit) where T <: AbstractFloat
    @unpack benchmark_subroutines = options
    @unpack total_steps = timer

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