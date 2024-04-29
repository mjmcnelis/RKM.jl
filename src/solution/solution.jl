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
    """Number of time steps taken"""
    time_steps_taken::MVector{1,Int64}
    """Number of function evaluations"""
    FE::MVector{1,Int64}
    """Number of Jacobian evaluations"""
    JE::MVector{1,Int64}
    """Step rejection rate (percentage)"""
    rejection_rate::MVector{1,Float64}
    """Runtime of ODE solver (excludes configuration) [seconds]"""
    runtime::MVector{1,Float64}
    """Memory used to store solution set (t,y) [bytes]"""
    solution_size::MVector{1,Int64}
    """Memory used to configure the solver [bytes]"""
    config_memory::MVector{1,Int64}
    """Excess memory used in time evolution loop [bytes]"""
    excess_memory::MVector{1,Int64}
    """Number of dynamical variables (assumed to be fixed)"""
    dimensions::MVector{1,Int64}
    """Float precision type"""
    precision::Type{T}
end

"""
    Solution(precision::Type{T} = Float64) where T <: AbstractFloat

Outer constructor for `Solution`.

Required parameters: `precision`
"""
function Solution(precision::Type{T} = Float64) where T <: AbstractFloat
    y = Vector{precision}()
    t = Vector{precision}()
    time_steps_taken = MVector{1,Int64}(0)
    FE = MVector{1,Int64}(0)
    JE = MVector{1,Int64}(0)
    rejection_rate = MVector{1,Float64}(0.0)
    runtime = MVector{1,Float64}(0.0)
    solution_size = MVector{1,Int64}(0)
    config_memory = MVector{1,Int64}(0)
    excess_memory = MVector{1,Int64}(0)
    dimensions = MVector{1,Int64}(0.0)

    return Solution(y, t, time_steps_taken, FE, JE, rejection_rate, runtime,
                    solution_size, config_memory, excess_memory, dimensions, precision)
end

function clear_solution!(sol::Solution)
    @unpack y, t, time_steps_taken, FE, JE, rejection_rate,
            runtime, solution_size, config_memory, excess_memory = sol
    empty!(y)
    empty!(t)
    sizehint!(y, 0)
    sizehint!(t, 0)
    time_steps_taken .= 0
    FE .= 0
    JE .= 0
    rejection_rate .= 0.0
    runtime .= 0.0
    solution_size .= 0
    config_memory .= 0
    excess_memory .= 0
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
    @unpack y, t, dimensions = sol
    # TODO: replace length(t) if use deleteat for PDEs
    y = reshape(y, dimensions[1], length(t)) |> transpose
    return y, t
end
