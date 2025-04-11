
function output_solution!(sol::Solution{T}, update_cache::UpdateCache{T},
                          options::SolverOptions{T}, t::Vector{T}) where T <: AbstractFloat
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