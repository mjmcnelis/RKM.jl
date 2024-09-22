
# TODO: replace ODEWrapper
struct ODEWrapperState{T, F} <: Wrapper where {T <: AbstractFloat, F <: Function}
    t::Vector{T}
    p::Vector{T}
    # make abstract_params::P that is separate from p
    dy_dt!::F
end

struct ODEWrapperParam{T, F} <: Wrapper where {T <: AbstractFloat, F <: Function}
    t::Vector{T}
    y::Vector{T}    # make a new vector or reuse one?
    dy_dt!::F
end

(ode_wrap!::ODEWrapperState)(f, y) = ode_wrap!.dy_dt!(f, y; p = ode_wrap!.p, t = ode_wrap!.t)
(ode_wrap!::ODEWrapperParam)(f, p) = ode_wrap!.dy_dt!(f, ode_wrap!.y; p, t = ode_wrap!.t)

# just make a rudimentary function that outputs something
# worry about organization, efficiency later
function post_sensitivity_analysis(sol, options, dy_dt!, p;
                                   # TODO: make SensitivityMethod struct
                                   jacobian_method = FiniteJacobian()
                                  )
    y, t = get_solution(sol)
    @unpack precision = options

    nt = length(t)
    ny = sol.dimensions[1]
    np = length(p)

    y0 = y[1,:]
    t0 = sol.t[1]

    y_tmp = copy(y0)

    # add to cache and generalize to matrices
    S = zeros(precision, ny)        # sensitivity cache
    S_tmp = zeros(precision, ny)

    Jy = zeros(precision, ny, ny)   # reuse/rename J in cache
    Jp = zeros(precision, ny, np)   # make new

    # FiniteDiff cache
    cache_y = JacobianCache(y0)
    cache_p = JacobianCache(y0)

    # create wrappers wrt y and p
    ode_wrap_y! = ODEWrapperState([t0], p, dy_dt!)
    ode_wrap_p! = ODEWrapperParam([t0], y0, dy_dt!)

    # sol.obj, sol.objS

    # ideally want to store in linear colummn format
    # and reshape it as a rank-3 tensor (Nt x Ny x Np)
    # TODO: call it sol.S
    yp = zeros(precision, ny)       # assume 1 parameter for now
    sizehint!(yp, nt*ny*np)         # minus -1?

    # looks like I can append a matrix to a vector (inner loop is rows)
    #=
    a = Float64[]
    b = [1.0 2.0; 3.0 4.0]
    append!(a, b')    # a = [1.0, 2.0, 3.0, 4.0]
    =#

    # Backward Euler
    for n in 1:nt-1
        y_tmp .= view(y, n+1, :)
        ode_wrap_p!.y .= y_tmp

        finite_difference_jacobian!(Jy, ode_wrap_y!, y_tmp, cache_y)
        finite_difference_jacobian!(Jp, ode_wrap_p!, p, cache_p)

        dt = t[n+1] - t[n]

        # TODO: add broadcast
        A = 1.0                 # stage coefficient (BackwardEuler)
        Jy .*= (-A*dt)          # J <- I - A.dt.J
        for i in diagind(Jy)
            Jy[i] += 1.0
        end

        S_tmp .= S .+ dt.*Jp        # option 1
        S_tmp .= inv(Jy)*S_tmp

        # S .+= dt.*Jp              # option 2
        # ldiv!(S_tmp, lu(Jy), S)   # need to add ldiv!

        # option 3: use LinearSolve

        append!(yp, S_tmp)
        S .= S_tmp
    end

    return yp
end