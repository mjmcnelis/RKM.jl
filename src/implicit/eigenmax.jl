abstract type EigenMaxMethod end

struct NoEigenMax <: EigenMaxMethod end

# use linear algebra eigvals, eigmax, only use for small systems
# TODO: how do I sort out real part
@kwdef struct LinearEigenMax <: EigenMaxMethod
    lambda_LR::Vector{ComplexF64} = ComplexF64[0.0]
end

# TODO: also store max eigenvalue (can you pass it to eigsolve?)
# should I store x_tmp in some manner (b/c I may want to reset x0 or update it)
# x_tmp = x0     # initialize before iteration loop
# x_tmp = x[idx] # update during iteration loop
# x0 <- x_tmp or x0 <- x0 # update after completed time step

@kwdef struct KrylovEigenMax <: EigenMaxMethod
    lambda_LR::Vector{ComplexF64} = ComplexF64[0.0]
    x0::Vector{ComplexF64} = ComplexF64[]   # maybe Vector{Vector{ComplexF64}}() is better?
    tol::Float64 = 1e-4 # note: varies problem to problem
    maxiter::Int64 = 10
    krylovdim::Int64 = 100 # check that <= ny
    verbosity::Int64 = 1
end

function reconstruct_eigenmax(eigenmax::NoEigenMax, args...)
    return eigenmax
end

function reconstruct_eigenmax(eigenmax::LinearEigenMax, args...)
    return eigenmax
end

function reconstruct_eigenmax(eigenmax::KrylovEigenMax, ny::Int64)
    # TODO: should imaginary part be random or zero?
    # @set! eigenmax.x0 = rand(ComplexF64, ny)
    @set! eigenmax.x0 = ComplexF64.(rand(ny))
    return eigenmax
end

# put methods here
function compute_max_eigenvalue!(::NoEigenMax, args...)
    return nothing
end

# note: t is tmp
function compute_max_eigenvalue!(eigenmax::LinearEigenMax, t::Vector{T},
             J::Union{Matrix{T}, SparseMatrixCSC{T,Int64}}) where T <: AbstractFloat

    @unpack lambda_LR = eigenmax

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

    # @show t lambda_LR[1]
    # println("")
    # sleep(0.1)
    # q()
    return nothing
end

function compute_max_eigenvalue!(eigenmax::KrylovEigenMax, t::Vector{T},
             J::Union{Matrix{T}, SparseMatrixCSC{T,Int64}}) where T <: AbstractFloat

    @unpack lambda_LR, x0, tol, maxiter, krylovdim, verbosity = eigenmax

    lambda, x, status = eigsolve(J, x0, 1, :LR; tol, maxiter, krylovdim, verbosity)
    idx = findfirst(x -> x == argmax(real, lambda), lambda)

    # if status.converged > 0
    #     @show t[1] lambda[idx] #eigmax(Matrix(J))
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