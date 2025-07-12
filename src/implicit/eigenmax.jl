abstract type EigenMaxMethod end

struct NoEigenMax <: EigenMaxMethod end

# note: only use for small systems
struct LinearEigenMax <: EigenMaxMethod end

# TODO: also store max eigenvalue (can you pass it to eigsolve?)
# should I store x_tmp in some manner (b/c I may want to reset x0 or update it)
# x_tmp = x0     # initialize before iteration loop
# x_tmp = x[idx] # update during iteration loop
# x0 <- x_tmp or x0 <- x0 # update after completed time step
@kwdef struct KrylovEigenMax <: EigenMaxMethod
    tol::Float64 = 1e-4 # note: varies problem to problem
    maxiter::Int64 = 10
    krylovdim::Int64 = 100 # check that <= ny
    verbosity::Int64 = 1
end

# put methods here
function compute_max_eigenvalue!(::NoEigenMax, args...)
    return nothing
end

# note: t is tmp
# TODO: pass update_cache and unpack lambda_LR, x0, J, y_tmp, f_tmp
#       (but may still need to pass (lambda_LR, x0) or (lambda_LR_tmp, x0_tmp), etc)
function compute_max_eigenvalue!(eigenmax::LinearEigenMax, lambda_LR::Vector{ComplexF64},
                                 x0::Vector{ComplexF64},
                                 J::Union{Matrix{T}, SparseMatrixCSC{T,Int64}},
                                 state_jacobian::JacobianMethod,
                                 ode_wrap!::ODEWrapperState, y::Vector{T},
                                 f::Vector{T}) where T <: AbstractFloat

    evaluate_jacobian!(state_jacobian, J, ode_wrap!, y, f)

    # compute eigenvalues
    if J isa SparseMatrixCSC
        @warn "LinearEigenMax() replaces sparse jacobian J with Matrix(J) for
               computing eigenvalues (only use for small systems)" maxlog = 1
        lambda = eigvals(Matrix(J))
    else
        lambda = eigvals(J)
    end

    # get eigenvalue with largest real part (still keeps imaginary part)
    lambda_LR[1] = argmax(real, lambda)

    # @show lambda_LR[1]
    # println("")
    # sleep(0.1)
    # q()
    return nothing
end

function compute_max_eigenvalue!(eigenmax::KrylovEigenMax, lambda_LR::Vector{ComplexF64},
                                 x0::Vector{ComplexF64},
                                 J::Union{Matrix{T}, SparseMatrixCSC{T,Int64}},
                                 state_jacobian::JacobianMethod,
                                 ode_wrap!::ODEWrapperState, y::Vector{T},
                                 f::Vector{T}) where T <: AbstractFloat

    evaluate_jacobian!(state_jacobian, J, ode_wrap!, y, f)

    tol = eigenmax.tol
    maxiter = eigenmax.maxiter
    krylovdim = eigenmax.krylovdim
    verbosity = eigenmax.verbosity

    lambda, x, status = eigsolve(J, x0, 1, :LR; tol, maxiter, krylovdim, verbosity)
    idx = argmax(real.(lambda))

    # if status.converged > 0
    #     @show t[1] lambda[idx]
    #     println("")
    #     sleep(0.1)
    #     # q()
    # end

    # get eigenvalue with largest real part (still keeps imaginary part)
    lambda_LR[1] = lambda[idx]

    # get corresponding eigenvector
    # TODO: use x_prev instead?
    x0 .= x[idx]

    return nothing
end