function compute_z0!(z0::SparseMatrixCSC{T,Int64},
                     A::SparseMatrixCSC{T,Int64}, t::T, t0::T) where T <: AbstractFloat
    @.. z0.nzval = A.nzval * (-(t - t0))
    return nothing
end

function compute_z0!(z0::Matrix{T}, A::Matrix{T}, t::T, t0::T) where T <: AbstractFloat
    @.. z0 = A * (-(t - t0))
    return nothing
end

# first prototype (assumes constant Jacobian)
function post_generator(sol::Solution{T}, options::SolverOptions{T}, dy_dt!::Function,
                        A::Union{Matrix{T}, SparseMatrixCSC{T,Int64}}, p::Vector{T};
                        abstract_params = nothing, dn::Int64 = 10, skip::Int64 = 100,
                        order::Int64 = 4) where T <: AbstractFloat

    @assert 0 <= order <= 4 "order = $order must be between 0 and 4"

    # TODO: evaluate and store Jacobian (call it A for now)
    precision = options.precision
    sensitivity = options.sensitivity

    t0 = sol.t[1]
    (; nt, ny, np) = get_dimensions(sol)

    # dn = stride used for finite difference
    # skip = interval between time points
    t_idxs = 2dn+1:skip:nt-2dn

    # Green's function
    # note: exponential function error if A is sparse

    A_dense = Matrix(A)
    @info "cond(A) = $(cond(A_dense))"

    G(t, t0) = exp(A_dense*(t - t0))      # or exp(-z)

    # allocate arrays
    y = zeros(precision, ny)
    f = zeros(precision, ny)
    z0 = copy(A)

    dydp = zeros(precision, ny, np)
    dfdp = zeros(precision, ny, np)
    # wonder if it's better to store subset of time points for dfdp
    # dfdp = zeros(precision, ny, np*5)
     # what about ny x np*order?
    dfdpdot = zeros(precision, ny, np, order)

    dSG = zeros(precision, ny, np)

    # should try to hug initial state better
    # SCE0 = zeros(precision, ny, np)

    # TODO: would be better to have configured already
    p = p .|> precision
    ode_wrap_p! = ODEWrapperParam([t0], y, abstract_params, dy_dt!)
    sensitivity = reconstruct_sensitivity(sensitivity, ode_wrap_p!, f, p, false)
    # param_jacobian should be ForwardJacobian (otherwise expect noise)
    param_jacobian = sensitivity.param_jacobian

    SG = Vector{precision}()
    sizehint!(SG, ny*np*length(t_idxs))

    initialized = false

    minusminus = [0.0, 0.0, -1.0, 1.0]  # n - 2dn
    minus = [-1.0, 1.0, 2.0, -4.0]      # n - dn
    curr = [0.0, -2.0, 0.0, 6.0]        # n
    plus = [1.0, 1.0, -2.0, -4.0]       # n + dn
    plusplus = [0.0, 0.0, 1.0, 1.0]     # n + 2dn

    denom_factor = [2.0, 1.0, 2.0, 1.0]

    for n in t_idxs
        t = sol.t[n]
        #=@time=# G0 = G(t, t0)
        #=@time=# compute_z0!(z0, A, t, t0)

        # can you reuse a matrix factorization if A is sparse?
        # note: don't use lu! unless reset A
        #=@time=# F = lu(A_dense)

        # assumes constant intervals
        dt = sol.t[n+dn] -  sol.t[n]

        # reset time derivatives (central differences)
        @.. dfdpdot = 0.0
        # dfdpdot1 = dfdp_p - dfdp_m
        # dfdpdot2 = dfdp_p - 2dfdp + dfdp_m
        # dfdpdot3 = dfdp_pp - 2dfdp_p + 2dfdp_m - dfdp_mm
        # dfdpdot4 = dfdp_pp - 4dfdp_p + 6dfdp - 4dfdp_m + dfdp_mm

        # evaluate dfdp at n-2dn
        t = sol.t[n-2dn]
        y .= view(sol.y, (1 + (n-1-2dn)*ny):(n-2dn)*ny)
        set_wrapper!(ode_wrap_p!, t, y)
        evaluate_jacobian!(param_jacobian, dfdp, ode_wrap_p!, p, f)
        # iterate 3rd and 4th derivatives
        for q in 1:order
            # continue if minusminuis[q] = 0
            @.. dfdpdot[:,:,q] += minusminus[q] * dfdp
        end

        # evaluate dfdp at n-dn
        t = sol.t[n-dn]
        y .= view(sol.y, (1 + (n-1-dn)*ny):(n-dn)*ny)
        set_wrapper!(ode_wrap_p!, t, y)
        evaluate_jacobian!(param_jacobian, dfdp, ode_wrap_p!, p, f)
        for q in 1:order
            @.. dfdpdot[:,:,q] += minus[q] * dfdp
        end

        # evaluate dfdp at n+dn
        t = sol.t[n+dn]
        y .= view(sol.y, (1 + (n-1+dn)*ny):(n+dn)*ny)
        set_wrapper!(ode_wrap_p!, t, y)
        evaluate_jacobian!(param_jacobian, dfdp, ode_wrap_p!, p, f)
        for q in 1:order
            @.. dfdpdot[:,:,q] += plus[q] * dfdp
        end

        # evaluate dfdp at n+2n
        t = sol.t[n+2dn]
        y .= view(sol.y, (1 + (n-1+2dn)*ny):(n+2dn)*ny)
        set_wrapper!(ode_wrap_p!, t, y)
        evaluate_jacobian!(param_jacobian, dfdp, ode_wrap_p!, p, f)
        for q in 1:order
            @.. dfdpdot[:,:,q] += plusplus[q] * dfdp
        end

        # current (n)
        # note: compute current dfdp last so can reuse in S0
        t = sol.t[n]
        y .= view(sol.y, (1 + (n-1)*ny):n*ny)
        set_wrapper!(ode_wrap_p!, t, y)
        evaluate_jacobian!(param_jacobian, dfdp, ode_wrap_p!, p, f)
        for q in 1:order
            @.. dfdpdot[:,:,q] += curr[q] * dfdp
        end

        # finite difference denominators
        for q in 1:order
            @.. dfdpdot[:,:,q] /= (denom_factor[q] * dt^q)
        end

        # note: S0 = dfdp, etc are just renaming of variables
        S0 = dfdp               # -A^{-1} dfdp
        @.. S0 *= -1.0
        ldiv!(F, S0)

        # non-allocating but takes a decent chunk of time
        # can do you them all at once initially?
        #=@time=# for q in 1:order
            dS = view(dfdpdot, :, :, q)    # -A^{-q+1} dfdpdotq
            @.. dS *= -1.0
            for _ in 1:q+1
                ldiv!(F, dS)
            end
        end

        # TODO: finish working out calcs for order < 4
        dS1 = view(dfdpdot, :, :, 1)
        dS2 = view(dfdpdot, :, :, 2)
        dS3 = view(dfdpdot, :, :, 3)
        dS4 = view(dfdpdot, :, :, 4)

        # is there a corresponding ODE for Gamma that I can project onto?
        # maybe its form is simpler to deal with

        # Gamma1 = I - G0
        # Gamma2 = I - G0*(I + z0)
        # Gamma3 = I - G0*(I + z0 + z0^2/2.0)
        # Gamma4 = I - G0*(I + z0 + z0^2/2.0 + z0^3/6.0)
        # Gamma5 = I - G0*(I + z0 + z0^2/2.0 + z0^3/6.0 + z0^4/24.0)

        # @.. dydp = 0.0
        # mul!(dSG, Gamma1, S0)
        # @.. dydp += dSG
        # mul!(dSG, Gamma2, dS1)
        # @.. dydp += dSG
        # mul!(dSG, Gamma3, dS2)
        # @.. dydp += dSG
        # mul!(dSG, Gamma4, dS3)
        # @.. dydp += dSG
        # mul!(dSG, Gamma5, dS4)
        # @.. dydp += dSG

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

        # note: integration by parts method not working (still a bug maybe?)
        # if !initialized
        #     @.. SCE0 = 0.0
        #     @.. SCE0 = S0 #+ dS1
        #     initialized = true
        # end
        # mul!(dSG, G0, SCE0)
        # @.. dydp = S0 #+ dS1
        # @.. dydp -= dSG

        append!(SG, dydp)
        # println("")
    end

    SG = reshape(SG, ny*np, length(t_idxs)) |> transpose

    return SG, t_idxs
end