# TODO: make docstring
function interpolate_solution(options::SolverOptions, sol::Solution,
                              precision::Type{T} = Float64;
                              dt_dense::Float64) where T <: AbstractFloat
    @unpack interpolator = options
    return interpolate_solution(interpolator, sol, precision; dt_dense)
end
