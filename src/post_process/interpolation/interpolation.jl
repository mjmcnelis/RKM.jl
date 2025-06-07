# TODO: make docstring
function interpolate_solution(options::SolverOptions, sol::Solution; dt_dense::Float64)
    @unpack interpolator, method, precision = options
    return interpolate_solution(interpolator, sol, method, precision; dt_dense)
end
