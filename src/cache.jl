
struct UpdateCache{T <: AbstractFloat}
    dy::Matrix{T}
    dy_LM::Matrix{T}
    y::Vector{T}
    y_tmp::Vector{T}
    f_tmp::Vector{T}
    f::Vector{T}
    J::Union{Matrix{T}, SparseMatrixCSC{T,Int64}}
    y1::Vector{T}
    y2::Vector{T}
    res::Vector{T}
    S::Matrix{T}
    S_tmp::Matrix{T}
    dS::Array{T,3}
    lambda_LR::Vector{ComplexF64}
    x0::Vector{ComplexF64}
    e_prev::Vector{T}
    tol_prev::Vector{T}
    dt_prev::Vector{T}
end

function UpdateCache(precision::Type{T}, y::Vector{T}, method::ODEMethod,
                     adaptive::AdaptiveTimeStep, dimensions::Int64, coefficients::Int64,
                     sensitivity::SensitivityMethod, state_jacobian::JacobianMethod,
                     eigenmax::EigenMaxMethod) where T <: AbstractFloat

    iteration = method.iteration
    stages = method.stages

    no_sensitivity = sensitivity isa NoSensitivity

    ny = dimensions                                             # state
    np = no_sensitivity ? 0 : coefficients                      # parameters
    nJ = (iteration isa Explicit && no_sensitivity) ? 0 : ny    # Jacobian
    nm = adaptive isa Fixed ? 0 : ny                            # primary/embedded
    # TODO: okay use res for both root solver and stepsize control?
    ne = iteration isa Explicit && adaptive isa Fixed ? 0 : ny  # residual error (res)
    nl = eigenmax isa NoEigenMax ? 0 : 1                        # eigenvalue

    if method isa LinearMultistep
        start_method = method.start_method
        dy = zeros(precision, ny, start_method.stages)
        dy_LM = zeros(precision, ny, stages)
    else
        dy = zeros(precision, ny, stages)
        dy_LM = Array{precision}(undef, 0, 0)
    end

    y_tmp = zeros(precision, ny)
    f_tmp = zeros(precision, ny)
    f = zeros(precision, ny)
    # TODO: may be better to split jacobian methods into sparse and non-sparse
    if hasproperty(state_jacobian, :sparsity) && size(state_jacobian.sparsity) == (ny, ny)
        J = state_jacobian.sparsity
        # TODO: still didn't work for Double64
        # J = state_jacobian.sparsity .|> precision
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
        # TODO: should imaginary part be random or zero?
        # x0 = rand(ComplexF64, ny)
        x0 = ComplexF64.(rand(ny))
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
