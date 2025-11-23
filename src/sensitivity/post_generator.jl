
# first prototype (assumes constant Jacobian)
function post_generator(sol::Solution{T}, options::SolverOptions{T},
                        dy_dt!::Function, A::Matrix{T}, p::Vector{T};
                        abstract_params = nothing) where T <: AbstractFloat

    # TODO: evaluate and store Jacobian (call it A for now)
    precision = options.precision
    sensitivity = options.sensitivity

    t0 = sol.t[1]
    (; nt, ny, np) = get_dimensions(sol)

    # tmp hard-code indices
    dn = 10                      # stride used for finite difference
    skip = 100                   # skipping time points is a big advantage
    t_idxs = 2dn+1:skip:nt-2dn

    # Green's function
    # note: exponential function error if A is sparse
    @info "cond(A) = $(cond(A))"
    G(t, t0) = exp(A*(t - t0))      # or exp(-z)
    z(t, t0) = -A*(t - t0)

    # allocate arrays
    y = zeros(precision, ny)
    f = zeros(precision, ny)

    dydp = zeros(precision, ny, np)
    dfdp = zeros(precision, ny, np)
    # wonder if it's better to store subset of time points for dfdp
    # dfdp = zeros(precision, ny, np*5)
    dSG = zeros(precision, ny, np)

    # TODO: include truncation order
    dfdpdot1 = zeros(precision, ny, np)
    dfdpdot2 = zeros(precision, ny, np)
    dfdpdot3 = zeros(precision, ny, np)
    dfdpdot4 = zeros(precision, ny, np)

    # TODO: would be better to have configured already
    p = p .|> precision
    ode_wrap_p! = ODEWrapperParam([t0], y, abstract_params, dy_dt!)
    sensitivity = reconstruct_sensitivity(sensitivity, ode_wrap_p!, f, p, false)
    # param_jacobian should be ForwardJacobian (otherwise expect noise)
    param_jacobian = sensitivity.param_jacobian

    SG = Vector{precision}()
    sizehint!(SG, ny*np*length(t_idxs))

    for n in t_idxs
        t = sol.t[n]
        #=@time=# G0 = G(t, t0)
        #=@time=# z0 = z(t, t0)

        # can you reuse a matrix factorization if A is sparse?
        #=@time=# F = lu(A)               # note: don't use lu! unless reset A

        # assumes constant intervals
        dt = sol.t[n+dn] -  sol.t[n]

        # reset time derivatives (central differences)
        @.. dfdpdot1 = 0.0          # dfdp_p - dfdp_m
        @.. dfdpdot2 = 0.0          # dfdp_p - 2dfdp + dfdp_m
        @.. dfdpdot3 = 0.0          # dfdp_pp - 2dfdp_p + 2dfdp_m - dfdp_mm
        @.. dfdpdot4 = 0.0          # dfdp_pp - 4dfdp_p + 6dfdp - 4dfdp_m + dfdp_mm

        # evaluate dfdp at n-2dn
        t = sol.t[n-2dn]
        y .= view(sol.y, (1 + (n-1-2dn)*ny):(n-2dn)*ny)
        set_wrapper!(ode_wrap_p!, t, y)
        evaluate_jacobian!(param_jacobian, dfdp, ode_wrap_p!, p, f)
        # iterate 3rd and 4th derivatives
        @.. dfdpdot3 -= dfdp
        @.. dfdpdot4 += dfdp

        # evaluate dfdp at n-dn
        t = sol.t[n-dn]
        y .= view(sol.y, (1 + (n-1-dn)*ny):(n-dn)*ny)
        set_wrapper!(ode_wrap_p!, t, y)
        evaluate_jacobian!(param_jacobian, dfdp, ode_wrap_p!, p, f)
        # iterate 1st - 4th derivatives
        @.. dfdpdot1 -= dfdp
        @.. dfdpdot2 += dfdp
        @.. dfdpdot3 += 2.0*dfdp
        @.. dfdpdot4 -= 4.0*dfdp

        # evaluate dfdp at n+dn
        t = sol.t[n+dn]
        y .= view(sol.y, (1 + (n-1+dn)*ny):(n+dn)*ny)
        set_wrapper!(ode_wrap_p!, t, y)
        evaluate_jacobian!(param_jacobian, dfdp, ode_wrap_p!, p, f)
        # iterate 1st - 4th derivatives
        @.. dfdpdot1 += dfdp
        @.. dfdpdot2 += dfdp
        @.. dfdpdot3 -= 2.0*dfdp
        @.. dfdpdot4 -= 4.0*dfdp

        # evaluate dfdp at n+2n
        t = sol.t[n+2dn]
        y .= view(sol.y, (1 + (n-1+2dn)*ny):(n+2dn)*ny)
        set_wrapper!(ode_wrap_p!, t, y)
        evaluate_jacobian!(param_jacobian, dfdp, ode_wrap_p!, p, f)
        # iterate 3rd and 4th derivatives
        @.. dfdpdot3 += dfdp
        @.. dfdpdot4 += dfdp

        # current (n)
        # note: compute current dfdp last so can reuse in S0
        t = sol.t[n]
        y .= view(sol.y, (1 + (n-1)*ny):n*ny)
        set_wrapper!(ode_wrap_p!, t, y)
        evaluate_jacobian!(param_jacobian, dfdp, ode_wrap_p!, p, f)
        # iterate 2nd and 4th time derivatives
        @.. dfdpdot2 -= 2.0*dfdp
        @.. dfdpdot4 += 6.0*dfdp

        # denominators
        @.. dfdpdot1 /= (2.0*dt)
        @.. dfdpdot2 /= dt^2
        @.. dfdpdot3 /= (2.0*dt^3)
        @.. dfdpdot4 /= dt^4

        # note: S0 = dfdp, etc are just renaming of variables
        S0 = dfdp               # -A^{-1} dfdp
        @.. S0 *= -1.0
        ldiv!(F, S0)

        dS1 = dfdpdot1          # -A^{-2} dfdpdot1
        @.. dS1 *= -1.0
        ldiv!(F, dS1)
        ldiv!(F, dS1)

        dS2 = dfdpdot2          # -A^{-3} dfdpdot2
        @.. dS2 *= -1.0
        ldiv!(F, dS2)
        ldiv!(F, dS2)
        ldiv!(F, dS2)

        dS3 = dfdpdot3          # -A^{-4} dfdpdot3
        @.. dS3 *= -1.0
        ldiv!(F, dS3)
        ldiv!(F, dS3)
        ldiv!(F, dS3)
        ldiv!(F, dS3)

        dS4 = dfdpdot4          # -A^{-5} dfdpdot4
        @.. dS4 *= -1.0
        ldiv!(F, dS4)
        ldiv!(F, dS4)
        ldiv!(F, dS4)
        ldiv!(F, dS4)
        ldiv!(F, dS4)

        # is there a corresponding ODE for Gamma that I can project onto?
        # maybe its form is simpler to deal with

        # Gamma1 = I - G0
        # Gamma2 = I - G0*(I + z0)
        # Gamma3 = I - G0*(I + z0 + z0^2/2.0)
        # Gamma4 = I - G0*(I + z0 + z0^2/2.0 + z0^3/6.0)
        # Gamma5 = I - G0*(I + z0 + z0^2/2.0 + z0^3/6.0 + z0^4/24.0)

        # dydp .= 0.0
        # mul!(dSG, Gamma1, S0)
        # dydp .+= dSG
        # mul!(dSG, Gamma2, dS1)
        # dydp .+= dSG
        # mul!(dSG, Gamma3, dS2)
        # dydp .+= dSG
        # mul!(dSG, Gamma4, dS3)
        # dydp .+= dSG
        # mul!(dSG, Gamma5, dS4)
        # dydp .+= dSG

        # not sure if I should keep going with this
        # z0dS1 = (t - t0) * (A \ ddtdfdp)
        # z0dS2 = (t - t0) * (A \ (A \ d2dt2dfdp))
        # z02dS2 = -(t - t0)^2 * (A \ d2dt2dfdp)
        # dydp0 = (I - G0) * (S0 + dS1 + dS2)
        # dydp1 = -G0*(z0dS1 + z0dS2 + 0.5*z02dS2)
        # dydp2 = 0*dfdp

        # dydp0 = -(I - G0) * (A \ (dfdp + A \ (ddtdfdp + A \ d2dt2dfdp)))
        # dydp1 = -(t - t0) * G0 * (A \ (ddtdfdp + A \ d2dt2dfdp))
        # dydp2 = 0.5 * (t - t0)^2 * G0 * (A \ d2dt2dfdp)

        # besides reducing projections, can you further
        # group dS terms to reduce linear solves?

        # current favorite b/c reduces projections
        @.. dydp = 0.0
        @.. dydp += S0 + dS1 + dS2 + dS3 + dS4

        #=@time=# dydp1 = -G0*(S0 + dS1 + dS2 + dS3 + dS4 +
                    z0*(dS1 + dS2 + dS3 + dS4 +
                        z0*((dS2 + dS3 + dS4)/2.0 +
                            z0*((dS3 + dS4)/6.0 +
                                z0*dS4/24.0
                                )
                            )
                        )
                    )
        @.. dydp += dydp1

        append!(SG, dydp)
        # println("")
    end

    SG = reshape(SG, ny*np, length(t_idxs)) |> transpose

    return SG, t_idxs
end