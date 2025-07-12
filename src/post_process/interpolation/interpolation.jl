# TODO: make docstring
function interpolate_solution(options::SolverOptions, sol::Solution; dt_dense::Float64)

    interpolator = options.interpolator
    method = options.method
    precision = options.precision

    return interpolate_solution(interpolator, sol, method, precision; dt_dense)
end
