abstract type EigenMaxMethod end

struct NoEigenMax <: EigenMaxMethod end

# use linear algebra eigvals, eigmax, only use for small systems
# TODO: how do I sort out real part
struct LinearEigenMax <: EigenMaxMethod end

# TODO: also store max eigenvalue (can you pass it to eigsolve?)
# should I store x_tmp in some manner (b/c I may want to reset x0 or update it)
# x_tmp = x0     # initialize before iteration loop
# x_tmp = x[idx] # update during iteration loop
# x0 <- x_tmp or x0 <- x0 # update after completed time step

@kwdef struct KrylovEigenMax <: EigenMaxMethod
    x0::Vector{ComplexF64} = ComplexF64[]   # maybe Vector{Vector{ComplexF64}}() is better?
    tol::Float64 = 1e-4 # note: varies problem to problem
    maxiter::Int64 = 10
    krylovdim::Int64 = 100 # check that <= ny
    # verbosity = 0
end

function reconstruct_eigenmax(eigenmax::NoEigenMax, args...)
    return eigenmax
end

function reconstruct_eigenmax(eigenmax::KrylovEigenMax, ny::Int64)
    # should imaginary part be random or zero?
    # @set! eigenmax.x0 = rand(ComplexF64, ny)
    @set! eigenmax.x0 = ComplexF64.(rand(ny))
    return eigenmax
end

# put methods here
function compute_max_eigenvalue!(::NoEigenMax, args...)
    return nothing
end

function compute_max_eigenvalue!(eigenmax::LinearEigenMax,
             J::Union{Matrix{T}, SparseMatrixCSC{T,Int64}}) where T <: AbstractFloat
    # TODO: fill out later
    if J isa SparseMatrixCSC
        lambda = 1
    else
        lambda = 1
        #
    end
    return nothing
end

function compute_max_eigenvalue!(eigenmax::KrylovEigenMax, t::Vector{T},
             J::Union{Matrix{T}, SparseMatrixCSC{T,Int64}}) where T <: AbstractFloat

    @unpack x0, tol, maxiter, krylovdim = eigenmax

    @time lambda, x, status = eigsolve(J, x0, 1, :LR; tol, maxiter, krylovdim)
    idx = findfirst(x -> x == argmax(real, lambda), lambda)
    @show t[1], lambda[idx]#, eigmax(Matrix(J))

    @show status
    if Bool(status.converged)
        sleep(0.1)
    end
    # TODO: use x_prev instead?
    x0 .= x[idx]
    # TODO: save eigenvalue
    return nothing
end