
# just make a playground function that outputs something
# worry about organization, efficiency later
function post_sensitivity_analysis(sol::Solution, options::SolverOptions,
                                   dy_dt!::Function, p::Vector{T};
                                   # TODO: make SensitivityMethod struct
                                   jacobian_method = FiniteJacobian()
                                  ) where {T <: AbstractFloat}
    t, y = get_solution(sol)
    @unpack precision = options

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
    # TODO: add sol.S (obj and objS)
    yp = zeros(precision, ny*np)
    sizehint!(yp, nt*ny*np)

    # Backward Euler
    for n in 1:nt-1
        # note: n+1 is specific to Backward Euler
        y_tmp .= view(y, n+1, :)            # can't do FastBroadcast here
        # set wrappers
        @.. ode_wrap_p!.y = y_tmp
        ode_wrap_y!.t[1] = t[n+1]
        ode_wrap_p!.t[1] = t[n+1]

        # compute Jacobian wrt y and p
        finite_difference_jacobian!(J, ode_wrap_y!, y_tmp, cache_y)
        finite_difference_jacobian!(S_tmp, ode_wrap_p!, p, cache_p)

        dt = t[n+1] - t[n]

        # TODO: add broadcast
        A = 1.0                 # stage coefficient (BackwardEuler)
        @.. J *= (-A*dt)        # J <- I - A.dt.J
        for k in diagind(J)
            J[k] = J[k] + 1.0
        end

        B = 1.0
        @.. S_tmp *= (B*dt)
        @.. S_tmp += S          # looks like stage calc order is reverse

        # any benefit in transposing the sensitivity ODE?
        # TODO: try using LinearSolve
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