"""
    nansafe_jacobian(y0::Vector{T}, t0::T1, dy_dt!::Function,
                     p::Vector{Float64} = Float64[];
                     chunk_size::Int64 = DEFAULT_CHUNK_THRESHOLD,
                     abstract_params = nothing) where {T <: AbstractFloat,
                                                       T1 <: AbstractFloat}

Estimates the sparsity pattern of the Jacobian J = df/dy via forward-mode
auto-differentiation. The Jacobian is evaluated at y = NaN, which only works if
NANSAFE_MODE_ENABLED = true in ForwardDiff.jl. Reducing the chunk size can help
remove excess matrix elements from the sparsity pattern. For more details
see https://juliadiff.org/ForwardDiff.jl/stable/user/advanced/.

Note: some nonzero matrix elements may be missed if `dy_dt!` removes traces of
      the dual number y = Dual(NaN, 1.0). For example, the following subroutines
      do this:

        z = isnan(y) ? zero(y) : y  # Dual(0.0, 0.0)
        z = max(y, 0.0)             # Dual(NaN, 0.0)
        z = min(0.0, y)             # Dual(NaN, 0.0)

Required parameters: `y0`, `t0`, `dy_dt!`

Optional parameters: `p`, `chunk_size`, `abstract_params`
"""
function nansafe_jacobian(y0::Vector{T}, t0::T1, dy_dt!::Function,
                          p::Vector{Float64} = Float64[];
                          chunk_size::Int64 = DEFAULT_CHUNK_THRESHOLD,
                          abstract_params = nothing) where {T <: AbstractFloat,
                                                            T1 <: AbstractFloat}
    print_nansafe_warning()

    ny = length(y0)
    J = zeros(ny, ny)

    y = [NaN for i in 1:ny]
    f = similar(y0)
    ode_wrap! = ODEWrapperState([t0], p, abstract_params, dy_dt!)
    cache = JacobianConfig(ode_wrap!, f, y, Chunk(chunk_size))

    jacobian!(J, ode_wrap!, f, y, cache)
    return sparse(J)
end

function nansafe_parameter_jacobian(y0::Vector{T}, t0::T1, dy_dt!::Function,
                                    p::Vector{Float64};
                                    chunk_size::Int64 = DEFAULT_CHUNK_THRESHOLD,
                                    abstract_params = nothing) where {T <: AbstractFloat,
                                                                      T1 <: AbstractFloat}
    print_nansafe_warning()

    ny = length(y0)
    np = length(p)
    Jp = zeros(ny, np)

    p_nan = [NaN for i in 1:np]

    # TODO: check if/where this matters
    # y = [NaN for i in 1:ny]
    y = copy(y0)
    f = similar(y0)

    ode_wrap_p! = ODEWrapperParam([t0], y, abstract_params, dy_dt!)
    cache = JacobianConfig(ode_wrap_p!, f, p_nan, Chunk(chunk_size))

    jacobian!(Jp, ode_wrap_p!, f, p_nan, cache)
    return sparse(Jp)
end

function print_nansafe_warning()
    if !NANSAFE_MODE_ENABLED
        @warn "nansafe mode is not enabled, sparsity pattern will likely be dense."
        println("""\nRun the following commands in your base (project) environment and
                restart the Julia REPL:

                using ForwardDiff, Preferences
                set_preferences!(ForwardDiff, "nansafe_mode" => true, force = true)

                The LocalPreferences.toml file can be edited directly in your base
                (project) environment (e.g. ~/.julia/environments/v1.10/).\n""")
    end
    return nothing
end