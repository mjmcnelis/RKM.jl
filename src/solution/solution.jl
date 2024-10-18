"""
Stores the solution vector `y(t)` of the ODE system in linear column format.

For example, the solution set `{y(0.0) = [1.0, 2.0, 3.0], y(0.5) = [4.0, 5.0, 6.0]}`
is stored as `y = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0], t = [0.0, 0.5]`.
"""
struct Solution{T <: AbstractFloat}
    """State vector of solution (stored as linear column)"""
    y::Vector{T}
    """First-order sensitivity coefficient matrix (stored as linear column)"""
    S::Vector{T}
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
    """Memory used to store sensitivity coefficients S [bytes]"""
    sensitivity_size::MVector{1,Int64}
    """Memory used to configure the solver [bytes]"""
    config_memory::MVector{1,Int64}
    """Excess memory used in time evolution loop [bytes]"""
    excess_memory::MVector{1,Int64}
    """Number of state variables (assumed to be fixed)"""
    dimensions::MVector{1,Int64}
    """Number of sensitivity coefficients/parameters"""
    coefficients::MVector{1,Int64}
end

"""
    Solution(precision::Type{T}) where T <: AbstractFloat

Outer constructor for `Solution`.

Required parameters: `precision`
"""
function Solution(precision::Type{T}) where T <: AbstractFloat
    y = Vector{precision}()
    S = Vector{precision}()
    t = Vector{precision}()
    time_steps_taken = MVector{1,Int64}(0)
    FE = MVector{1,Int64}(0)
    JE = MVector{1,Int64}(0)
    rejection_rate = MVector{1,Float64}(0.0)
    runtime = MVector{1,Float64}(0.0)
    sensitivity_size = MVector{1,Int64}(0)
    solution_size = MVector{1,Int64}(0)
    config_memory = MVector{1,Int64}(0)
    excess_memory = MVector{1,Int64}(0)
    dimensions = MVector{1,Int64}(0)
    coefficients = MVector{1,Int64}(0)

    # never understood why do {precision}
    return Solution{precision}(y, S, t, time_steps_taken, FE, JE, rejection_rate, runtime,
                               solution_size, sensitivity_size, config_memory,
                               excess_memory, dimensions, coefficients)
end

function clear_solution!(sol::Solution)
    @unpack y, S, t, time_steps_taken, FE, JE, rejection_rate,
            runtime, solution_size, sensitivity_size, config_memory,
            excess_memory, dimensions, coefficients = sol
    empty!(y)
    empty!(S)
    empty!(t)
    sizehint!(y, 0)
    sizehint!(S, 0)
    sizehint!(t, 0)
    time_steps_taken .= 0
    FE .= 0
    JE .= 0
    rejection_rate .= 0.0
    runtime .= 0.0
    solution_size .= 0
    sensitivity_size .= 0
    config_memory .= 0
    excess_memory .= 0
    dimensions .= 0
    coefficients .= 0
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
    return t, y
end

function get_sensitivity(sol::Solution)
    @unpack S, t, dimensions, coefficients = sol

    ny = dimensions[1]
    np = coefficients[1]
    # TODO: replace length(t) if use deleteat for PDEs
    S = reshape(S, ny*np, length(t)) |> transpose
    return t, S
end
