
# struct UpdateCache{T <: AbstractFloat, TE <: Union{T, Complex{T}}}

# const AFCAF{T} = Union{T, Complex{T}} where T <: AbstractFloat
# struct UpdateCache{T <: AbstractFloat, TE <: AFCAF{T}}

struct UpdateCache{T <: AbstractFloat, TCT <: Union{T, Complex{T}}}
    dy::Matrix{TCT}
    dy_LM::Matrix{T}
    y::Vector{TCT}
    y_tmp::Vector{TCT}
    f_tmp::Vector{TCT}
    f::Vector{TCT}
    J::Union{Matrix{T}, SparseMatrixCSC{T,Int64}}
    y1::Vector{T}
    y2::Vector{T}
    res::Vector{T}
    S::Matrix{T}
    S_tmp::Matrix{T}
    dS::Array{T,3}
    lambda_LR::Vector{ComplexF64}   # TODO: use Complex{T}
    x0::Vector{ComplexF64}
    e_prev::Vector{T}
    tol_prev::Vector{T}
    dt_prev::Vector{T}
end

function UpdateCache(precision::Type{T}, y::Vector{TCT}, method::ODEMethod,
                     adaptive::AdaptiveTimeStep, dimensions::Int64, coefficients::Int64,
                     sensitivity::SensitivityMethod, state_jacobian::JacobianMethod,
                     eigenmax::EigenMaxMethod) where {T <: AbstractFloat,
                                                      TCT <: Union{T, Complex{T}}}

    iteration = method.iteration
    stages = method.stages

    sparsity = state_jacobian.sparsity

    no_sensitivity = sensitivity isa NoSensitivity

    ny = dimensions                                             # state
    np = no_sensitivity ? 0 : coefficients                      # parameters
    nJ = (iteration isa Explicit && no_sensitivity) ? 0 : ny    # Jacobian
    nm = adaptive isa Fixed ? 0 : ny                            # primary/embedded
    ne = iteration isa Explicit && adaptive isa Fixed ? 0 : ny  # residual error (res)
    nl = eigenmax isa NoEigenMax ? 0 : 1                        # eigenvalue

    if method isa LinearMultistep
        start_method = method.start_method
        dy = zeros(TCT, ny, start_method.stages)
        dy_LM = zeros(TCT, ny, stages)
    else
        dy = zeros(TCT, ny, stages)
        dy_LM = Array{precision}(undef, 0, 0)
    end

    y_tmp = zeros(TCT, ny)
    f_tmp = zeros(TCT, ny)
    f = zeros(TCT, ny)

    if size(sparsity) == (ny, ny)
        J = sparsity .|> precision
    else
        J = zeros(precision, nJ, nJ)
    end

    y1 = zeros(precision, nm)
    y2 = zeros(precision, nm)
    res = zeros(precision, ne)

    S = zeros(precision, ny, np)
    S_tmp = zeros(precision, ny, np)
    dS = zeros(precision, ny, np, stages)

    lambda_LR = zeros(ComplexF64, nl)

    if eigenmax isa KrylovEigenMax
        random_imaginary = eigenmax.random_imaginary
        if random_imaginary
            x0 = rand(ComplexF64, ny)
        else
            x0 = ComplexF64.(rand(ny))
        end
    else
        x0 = ComplexF64[]
    end

    # for time step controller
    e_prev = ones(precision, 2)
    tol_prev = ones(precision, 2)
    dt_prev  = ones(precision, 3)

    return UpdateCache(dy, dy_LM, y, y_tmp, f_tmp, f, J, y1, y2, res, S,
                       S_tmp, dS, lambda_LR, x0, e_prev, tol_prev, dt_prev)
end
