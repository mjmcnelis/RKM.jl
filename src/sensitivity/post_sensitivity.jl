
# just make a playground function that outputs something
function post_sensitivity_analysis(sol::Solution, options::SolverOptions, dy_dt!::Function,
                                   p::Vector{T}) where {T <: AbstractFloat}

    t, y = get_solution(sol)

    precision = options.precision

    nt = length(t)
    ny = sol.dimensions[1]
    np = length(p)

    y0 = y[1,:]
    t0 = sol.t[1]

    y_tmp = copy(y0)

    S = zeros(precision, ny, np)        # add sensitivity cache
    S_tmp = zeros(precision, ny, np)

    J = zeros(precision, ny, ny)        # reuse J in cache

    # FiniteDiff cache
    cache_y = JacobianCache(y0)
    cache_p = JacobianCache(p, y0)

    # create wrappers wrt y and p
    abstract_params = nothing   # TMP
    ode_wrap_y! = ODEWrapperState([t0], p, abstract_params, dy_dt!)
    ode_wrap_p! = ODEWrapperParam([t0], y0, abstract_params, dy_dt!)

    # store as linear column, then reshape as nt x (ny*np) matrix
    yp = zeros(precision, ny*np)
    sizehint!(yp, nt*ny*np)

    # Backward Euler
    for n in 1:nt-1
        # note: n+1 is specific to Backward Euler
        @.. y_tmp = y[n+1,:]
        # set wrappers
        set_wrapper!(ode_wrap_y!, t[n+1])
        set_wrapper!(ode_wrap_p!, t[n+1], y_tmp)

        # compute Jacobian wrt y and p
        finite_difference_jacobian!(J, ode_wrap_y!, y_tmp, cache_y)
        finite_difference_jacobian!(S_tmp, ode_wrap_p!, p, cache_p)

        dt = t[n+1] - t[n]

        A = 1.0                 # stage coefficient (BackwardEuler)
        @.. J *= (-A*dt)        # J <- I - A.dt.J
        for k in diagind(J)
            J[k] = J[k] + 1.0
        end

        B = 1.0
        @.. S_tmp *= (B*dt)
        @.. S_tmp += S          # looks like stage calc order is reverse

        # any benefit in transposing the sensitivity ODE?
        F = lu!(J)
        ldiv!(F, S_tmp)

        # TMP for debugging reshape
        # S_tmp[1,2] = 0.01
        # S_tmp[2,1] = -0.01

        append!(yp, S_tmp)

        # comment if done debugging reshape
        # S_tmp[1,2] = 0.0
        # S_tmp[2,1] = 0.0

        @.. S = S_tmp
    end

    return yp
end

function psa_green_function(sol::Solution, options::SolverOptions,
                            dy_dt!::Function, p::Vector{T};
                            abstract_params = nothing) where {T <: AbstractFloat}

    t, y = get_solution(sol)

    precision = options.precision

    nt = length(t)
    ny = sol.dimensions[1]

    y0 = y[1,:]
    t0 = sol.t[1]

    y_tmp = copy(y0)
    z0_tmp = zeros(precision, ny, ny)
    J = zeros(precision, ny, ny)
    cache_y = JacobianCache(y0)
    ode_wrap_y! = ODEWrapperState([t0], p, abstract_params, dy_dt!)

    z0 = zeros(precision, ny^2)
    sizehint!(z0, nt*ny^2)

    # Backward Euler
    for n in 1:nt-1
        # note: n+1 is specific to Backward Euler
        @.. y_tmp = y[n+1,:]
        ode_wrap_y!.t[1] = t[n+1]

        # compute Jacobian
        finite_difference_jacobian!(J, ode_wrap_y!, y_tmp, cache_y)
        dt = t[n+1] - t[n]

        @.. z0_tmp += J*dt
        append!(z0, z0_tmp)
    end
    z0 = reshape(z0, ny^2, nt) |> transpose

    return z0
end