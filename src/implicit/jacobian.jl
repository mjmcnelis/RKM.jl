
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
