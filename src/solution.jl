"""
Stores the solution vector `y(t)` of the ODE system in linear column format.

For example, the solution set `{y(0.0) = [1.0, 2.0, 3.0], y(0.5) = [4.0, 5.0, 6.0]}`
is stored as `y = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0], t = [0.0, 0.5]`
"""
struct Solution{T <: AbstractFloat}
    """Solution vector (stored as linear column)"""
    y::Vector{T}
    """Time vector"""
    t::Vector{T}
    """Number of function evaluations"""
    FE::MVector{1,Int64}
    """Number of dynamical variables"""
    dimensions::Int64
end

"""
Outer constructor for `Solution`.

Required parameters: `precision`, `dimensions`.
"""
function Solution(; precision::Type{T}, dimensions::Int64) where T <: AbstractFloat
    y = Vector{precision}()
    t = Vector{precision}()
    FE = MVector{1,Int64}(0)

    Solution(y, t, FE, dimensions)
end

"""
    update_solution!(sol::Solution, y::Vector{T},
                     t::MVector{1,T2}) where {T <: AbstractFloat, T2 <: AbstractFloat}

Appends the state vector `y` at the current time `t` to the solution `sol`.

Required parameters: `sol`, `y`, `t`

Note: currently `y` and `t` can be different float types.
"""
function update_solution!(sol::Solution, y::Vector{T},
                          t::MVector{1,T2}) where {T <: AbstractFloat, T2 <: AbstractFloat}
    append!(sol.y, y)
    append!(sol.t, t[1])
    nothing
end

"""
    sizehint_solution!(sol::Solution, t_span::TimeSpan, dimensions::Int64)

Applies `sizehint!` to vector fields in the solution `sol`.

Required parameters: `sol`, `t_span`, `dimensions`
"""
function sizehint_solution!(sol::Solution, t_span::TimeSpan, dimensions::Int64)
    # TODO: whether or not I call this depends if save at regular intervals
    @unpack t0, tf, dt0 = t_span
    steps = round((tf - t0)/dt0) |> Int64
    # TODO: why is steps+1 not sufficient?        
    sizehint!(sol.y, dimensions*(steps + 2))
    sizehint!(sol.t, steps + 2)
    nothing
end

"""
    reshape_solution(sol::Solution)

Reshapes a solution vector `sol.y` of length `D*N` into a `N x D` matrix format 
(`N` = time steps, `D` = dimensions).

For example, the solution set `{y(0.0) = [1.0, 2.0, 3.0], y(0.5) = [4.0, 5.0, 6.0]}`
is stored in linear column format as `y = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]`.

The solution set is then reshaped as `y = [1.0 2.0 3.0; 4.0 5.0 6.0]`.
"""
function reshape_solution(sol::Solution)
    y, t = sol.y, sol.t
    # TODO: replace length(t) if use deleteat for PDEs
    y = reshape(y, sol.dimensions, length(t))'
    y, t
end
