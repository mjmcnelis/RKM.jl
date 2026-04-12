function compute_z0!(z0::SparseMatrixCSC{T,Int64},
                     A::SparseMatrixCSC{T,Int64}, t::T, t0::T) where T <: AbstractFloat
    @.. z0.nzval = A.nzval * (-(t - t0))
    return nothing
end

function compute_z0!(z0::Matrix{T}, A::Matrix{T}, t::T, t0::T) where T <: AbstractFloat
    @.. z0 = A * (-(t - t0))
    return nothing
end

# TODO: add struct options to use ODE solver or Krylov (so can take out kwargs)
function green_matrix_product!(GX::Matrix{T}, A::Union{Matrix{T}, SparseMatrixCSC{T,Int64}},
                               X::Matrix{T}, t::T, t0::T; verbosity::Int64 = 1,
                               krylovdim::Int64 = 10, tol::Float64 = 1e-12,
                               maxiter::Int64 = 100) where T <: AbstractFloat
    # X = ny x np matrix (constant)
    # GX = ny x np matrix
    # Gx = ny x 1 vector

    # note: KrylovKit doesn't support matrix-matrix product exp(A*(t-t0))*X
    ny, np = size(X)
    for j in 1:np
        x = view(X, :, j)
        # this is a lot slower than I thought
        Gx, info = exponentiate(A, t - t0, x; verbosity, krylovdim, tol, maxiter)
        # Gx = expv(t - t0, A, x; m = krylovdim)
        @show info
        @.. GX[:,j] = Gx
    end
    return nothing
end

# first prototype (assumes constant Jacobian)
function post_generator(sol::Solution{T}, options::SolverOptions{T}, dy_dt!::Function,
                        A::Union{Matrix{T}, SparseMatrixCSC{T,Int64}}, p::Vector{T};
                        abstract_params = nothing, dn::Int64 = 10, skip::Int64 = 100,
                        order::Int64 = 2) where T <: AbstractFloat

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
    # A_dense = Matrix(A)
    # @info "cond(A) = $(cond(A_dense))"
    # G(t, t0) = exp(A_dense*(t - t0))      # or exp(-z)

    # allocate arrays
    y = zeros(precision, ny)
    f = zeros(precision, ny)

    GX = zeros(precision, ny, np)
    z0 = copy(A)

    dydp = zeros(precision, ny, np)
    dfdp = zeros(precision, ny, np)
    # wonder if it's better to store subset of time points for dfdp
    # dfdp = zeros(precision, ny, np*5)
     # what about ny x np*order?
    dfdpdot = zeros(precision, ny, np, order)

    S_tmp = zeros(precision, ny, np)
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

    F = lu(A)

    t_prev = t0

    for n in t_idxs
        t = sol.t[n]
        compute_z0!(z0, A, t, t_prev)

        # can you reuse a matrix factorization if A is sparse?
        # note: don't use lu! unless reset A
        # #=@time=# if A isa SparseMatrixCSC
        #     lu!(F, A)
        # else
        #     F = lu(A)
        # end
        @btime lu!($F, $A)

        # assumes constant intervals
        dt = sol.t[n+dn] - sol.t[n]

    # @time begin
        # reset time derivatives (central differences)
        @.. dfdpdot = 0.0

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
    # end

    # @time begin
        # note: S0 = dfdp, etc are just renaming of variables
        S0 = dfdp               # -A^{-1} dfdp
        @.. S0 *= -1.0
        ldiv!(S_tmp, F, S0)
        @.. S0 = S_tmp

        # non-allocating but takes a decent chunk of time
        # can do you them all at once initially?
        for q in 1:order
            dS = view(dfdpdot, :, :, q)    # -A^{-q+1} dfdpdotq
            @.. dS *= -1.0
            for _ in 1:q+1
                ldiv!(S_tmp, F, dS)
                @.. dS = S_tmp
            end
        end

        # TODO: finish working out calcs for order < 4
        dS1 = view(dfdpdot, :, :, 1)
        dS2 = view(dfdpdot, :, :, 2)
    # end

        # it's more to do with the local error of the expansion as opposed to global
        # note: this isn't optimized but it takes up a small fraction
        #=@time=# dydp_tmp = -(S0 + dS1 + dS2 + z0*(dS1 + dS2 + z0*dS2/2.0))
        @.. dydp += dydp_tmp

        ##=@time=# green_matrix_product!(GX, A, dydp, t, t_prev)
        @btime green_matrix_product!($GX, $A, $dydp, $t, $t_prev)

        # reset dydp to GX
        @.. dydp = GX
        @.. dydp += S0 + dS1 + dS2

        append!(SG, dydp)
        # println("")

        t_prev = t

        break
    end

    # SG = reshape(SG, ny*np, length(t_idxs)) |> transpose

    return SG, t_idxs
end