
# TODO: should store max eigenvalue time series if computing it
function get_eigenvalues(sol::Solution{T}, dy_dt!::Function, options::SolverOptions{T},
                         p::Vector{Float64} = Float64[]; abstract_params = nothing,
                         dt_dense::Float64) where T <: AbstractFloat

    @info "Computing eigenvalues of the state jacobian"
    # TODO: prefer suppressing interpolation @info here
    t, y = interpolate_solution(options, sol; dt_dense)

    precision = options.precision

    nt = length(t)
    ny = sol.dimensions[1]

    t0 = t[1]
    y0 = copy(y[1,:])

    y_tmp = copy(y0)
    J = zeros(precision, ny, ny)
    cache = JacobianCache(y0)

    ode_wrap! = ODEWrapperState([t0], p, abstract_params, dy_dt!)

    lambda = Vector{ComplexF64}()
    sizehint!(lambda, nt*ny)

    for n in 1:nt
        @.. y_tmp = y[n,:]

        # assume finite diff for now
        finite_difference_jacobian!(J, ode_wrap!, y_tmp, cache)

        append!(lambda, eigvals(J))
    end

    lambda = reshape(lambda, ny, nt) |> transpose

    return t, lambda
end