
abstract type RootMethod end
struct Newton <: RootMethod end         # TODO: add line search later
struct FixedPoint <: RootMethod end

abstract type JacobianMethod end

@kwdef struct ForwardJacobian{JC} <: JacobianMethod where JC <: JacobianConfig
    cache::JC = JacobianConfig(nothing, [0.0], [0.0])
    evaluations::MVector{1,Int64} = MVector{1,Int64}(0)
end
@kwdef struct ForwardColorJacobian{JC, T} <: JacobianMethod where {JC <: ForwardColorJacCache,
                                                                   T <: AbstractFloat}
    cache::JC = ForwardColorJacCache(nothing, [0.0])
    sparsity::SparseMatrixCSC{T,Int64} = SparseMatrixCSC(Float64[;;])
    evaluations::MVector{1,Int64} = MVector{1,Int64}(0)
end

@kwdef struct FiniteJacobian{JC, T} <: JacobianMethod where {JC <: JacobianCache,
                                                             T <: AbstractFloat}
    cache::JC = JacobianCache([0.0])
    sparsity::SparseMatrixCSC{T,Int64} = SparseMatrixCSC(Float64[;;])
    evaluations::MVector{1,Int64} = MVector{1,Int64}(0)
end

abstract type StageFinder end

# good enough start (wrap caches later)
@kwdef struct ImplicitStageFinder{JM, AF} <: StageFinder where {JM <: JacobianMethod,
                                                                AF <: AbstractFactorization
                                                               }
    root_method::RootMethod = Newton()
    jacobian_method::JM     = FiniteJacobian()
    linear_method::AF       = LUFactorization()
     # TODO: reuse adaptive epsilon or 100x smaller?
    epsilon::Float64        = 1e-8
    max_iterations::Int64   = 10
    # TODO: make outer constructor to check p_norm value
    p_norm::Float64         = 2.0
    # add iterations_per_stage
end

function set_jacobian_cache(stage_finder::ImplicitStageFinder, ode_wrap!, f, y)
    # reset jacobian evaluations
    stage_finder.jacobian_method.evaluations .= 0

    # TODO: only set cache is use Newton method
    @unpack jacobian_method = stage_finder
    if jacobian_method isa FiniteJacobian
        @unpack sparsity = jacobian_method
        if all(size(sparsity) .== length(y))
            colorvec = matrix_colors(sparsity)
            cache = JacobianCache(y; colorvec, sparsity)
        else
            cache = JacobianCache(y)
        end
    elseif jacobian_method isa ForwardJacobian
        cache = JacobianConfig(ode_wrap!, f, y)
    elseif jacobian_method isa ForwardColorJacobian
        @unpack sparsity = jacobian_method
        if all(size(sparsity) .== length(y))
            colorvec = matrix_colors(sparsity)
            cache = ForwardColorJacCache(ode_wrap!, y; colorvec, sparsity)
        else
            cache = ForwardColorJacCache(ode_wrap!, y)
        end
    end
    @set! stage_finder.jacobian_method.cache = cache
    return stage_finder
end

function evaluate_system_jacobian!(jacobian_method::ForwardJacobian,
                                   J, ode_wrap!, y, f)
    @unpack cache, evaluations = jacobian_method
    jacobian!(J, ode_wrap!, f, y, cache)
    evaluations[1] += 1
    return nothing
end

function evaluate_system_jacobian!(jacobian_method::ForwardColorJacobian,
                                   J, ode_wrap!, y, args...)
    @unpack cache, evaluations = jacobian_method
    forwarddiff_color_jacobian!(J, ode_wrap!, y, cache)
    evaluations[1] += 1
    return nothing
end

function evaluate_system_jacobian!(jacobian_method::FiniteJacobian,
                                   J, ode_wrap!, y, args...)
    @unpack cache, evaluations = jacobian_method
    finite_difference_jacobian!(J, ode_wrap!, y, cache)
    evaluations[1] += 1
    return nothing
end

function nansafe_jacobian(y0::Vector{T}, t0::T1, dy_dt!::Function,
                          p::Vector{Float64} = Float64[];
                          abstract_params = nothing) where {T <: AbstractFloat,
                                                            T1 <: AbstractFloat}
    if !NANSAFE_MODE_ENABLED
        @warn "nansafe_mode is false, sparsity pattern will likely be dense."
        println("""\nRun the following commands in your base (project) environment and
                restart the Julia REPL:

                    using ForwardDiff, Preferences
                    set_preferences!(ForwardDiff, "nansafe_mode" => true, force = true)

                The LocalPreferences.toml file can be edited directly in your base
                (project) environment (e.g. ~/.julia/environments/v1.10/).\n""")
    end
    ny = length(y0)
    J = zeros(ny, ny)

    y = NaN.*y0
    f = similar(y0)
    ode_wrap! = ODEWrapperState([t0], p, abstract_params, dy_dt!)

    jacobian!(J, ode_wrap!, f, y)
    return sparse(J)
end
