"""
Stores the solution vector `y(t)` of the ODE system in linear column format.

For example, the solution set `{y(0.0) = [1.0, 2.0, 3.0], y(0.5) = [4.0, 5.0, 6.0]}`
is stored as `y = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0], t = [0.0, 0.5]`.
"""
struct Solution{T <: AbstractFloat}
    """State vector of solution (stored as linear column)"""
    y::Vector{T}
    """Time vector of solution"""
    t::Vector{T}
    """Number of function evaluations"""
    FE::MVector{1,Int64}
    """Number of Jacobian evaluations"""
    JE::MVector{1,Int64}
    """Step rejection rate (percentage)"""
    rejection_rate::MVector{1,Float64}
    """Runtime of ODE solver (excludes configuration)"""
    runtime::MVector{1,Float64}
    """Memory used to store solution set (t,y)"""
    memory_storage::Vector{String}
    """Excess memory used in time evolution loop"""
    excess_memory::Vector{String}
    """Excess number of allocations in time evolution loop"""
    excess_allocations::MVector{1,Int64}
    """Number of dynamical variables"""
    dimensions::Int64
end

"""
    Solution(; precision::Type{T}, dimensions::Int64) where T <: AbstractFloat

Outer constructor for `Solution`.

Required parameters: `precision`, `dimensions`
"""
function Solution(; precision::Type{T}, dimensions::Int64) where T <: AbstractFloat
    y = Vector{precision}()
    t = Vector{precision}()
    FE = MVector{1,Int64}(0)
    JE = MVector{1,Int64}(0)
    rejection_rate = MVector{1,Float64}(0.0)
    runtime = MVector{1,Float64}(0.0)
    memory_storage = [""]
    excess_memory = [""]
    excess_allocations = MVector{1,Int64}(0)

    return Solution(y, t, FE, JE, rejection_rate, runtime, memory_storage,
                    excess_memory, excess_allocations, dimensions)
end

"""
    sizehint_solution!(sol::Solution, t_range::TimeRange, dimensions::Int64)

Applies `sizehint!` to the vector fields `y` and `t` in the solution `sol`.

Required parameters: `sol`, `t_range`, `dimensions`
"""
function sizehint_solution!(sol::Solution, t_range::TimeRange, dimensions::Int64)
    # TODO: whether or not I call this depends if save at regular intervals
    @unpack t0, tf, dt0 = t_range
    steps = round((tf - t0)/dt0) |> Int64
    # TODO: why is steps+1 not sufficient?
    sizehint!(sol.y, dimensions*(steps + 2))
    sizehint!(sol.t, steps + 2)
    return nothing
end

"""
    get_solution(sol::Solution)

Returns the solution tuple `(y,t)` from `sol`. The solution vector `y`, which has a length
`D*N`, is reshaped into an `N x D` matrix (`N` = time steps, `D` = dimensions).

For example, the solution set `{y(0.0) = [1.0, 2.0, 3.0], y(0.5) = [4.0, 5.0, 6.0]}`
is stored in linear column format as `y = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]`. The solution
vector is then reshaped as `y = [1.0 2.0 3.0; 4.0 5.0 6.0]`.
"""
function get_solution(sol::Solution)
    y, t = sol.y, sol.t
    # TODO: replace length(t) if use deleteat for PDEs
    y = reshape(y, sol.dimensions, length(t))'
    return y, t
end

function get_stats(sol::Solution)
    println("time steps           = $(length(sol.t))")
    println("step rejection rate  = $(sol.rejection_rate[1]) %")
    println("function evaluations = $(sol.FE[1])")
    println("jacobian evaluations = $(sol.JE[1])")
    println("solver runtime       = $(round(sol.runtime[1], sigdigits = 4)) seconds")
    println("solution storage     = $(sol.memory_storage[1])")
    println("excess memory        = $(sol.excess_memory[1])")
    println("excess allocations   = $(sol.excess_allocations[1])")
end
