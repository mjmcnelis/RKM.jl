"""
Stores the solution vector `y(t)` of the ODE system in linear column format.

For example, the solution set `{y(0.0) = [1.0, 2.0, 3.0], y(0.5) = [4.0, 5.0, 6.0]}`
is stored as `y = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0], t = [0.0, 0.5]`.
"""
struct Solution{T <: AbstractFloat}
    """Time series of solution"""
    t::Vector{T}
    """State variables of solution (stored as linear column)"""
    y::Vector{T}
    """Time derivatives of solution (stored as linear column)"""
    f::Vector{T}
    """Intermediate stages of solution (stored as linear column)"""
    dy::Vector{T}
    """First-order sensitivity coefficients (stored as linear column)"""
    S::Vector{T}
    """Eigenvalue of state Jacobian with largest real part"""
    lambda_LR::Vector{ComplexF64}
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
    t = Vector{precision}()
    y = Vector{precision}()
    f = Vector{precision}()
    dy = Vector{precision}()
    S = Vector{precision}()
    lambda_LR = Vector{ComplexF64}()
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
    return Solution{precision}(t, y, f, dy, S, lambda_LR, time_steps_taken, FE, JE,
                               rejection_rate, runtime, solution_size, sensitivity_size,
                               config_memory, excess_memory, dimensions, coefficients)
end

function clear_solution!(sol::Solution)
    @unpack t, y, f, dy, S, lambda_LR, time_steps_taken, FE, JE, rejection_rate,
            runtime, solution_size, sensitivity_size, config_memory,
            excess_memory, dimensions, coefficients = sol
    empty!(t)
    empty!(y)
    empty!(f)
    empty!(dy)
    empty!(S)
    empty!(lambda_LR)
    sizehint!(t, 0)
    sizehint!(y, 0)
    sizehint!(f, 0)
    sizehint!(dy, 0)
    sizehint!(S, 0)
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

Returns the solution tuple `(t,y)` from `sol`. The solution vector `y`, which has a length
`D*N`, is reshaped into an `N x D` matrix (`N` = time steps, `D` = dimensions).

For example, the solution set `{y(0.0) = [1.0, 2.0, 3.0], y(0.5) = [4.0, 5.0, 6.0]}`
is stored in linear column format as `y = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]`. The solution
vector is then reshaped as `y = [1.0 2.0 3.0; 4.0 5.0 6.0]`.
"""
function get_solution(sol::Solution)
    @unpack t, y, dimensions = sol
    if isempty(y)
        error("State variables y = $y are empty, set save_solution = true")
    end
    ny = dimensions[1]
    nt = length(t)
    # TODO: replace nt if use deleteat for PDEs
    y_reshape = reshape(y, ny, nt) |> transpose
    return t, y_reshape
end

function get_time_derivative(sol::Solution)
    @unpack t, f, dimensions = sol
    if isempty(f)
        error("Time derivatives f = $f are empty, set save_solution = true \
               and save_time_derivative = true")
    end
    ny = dimensions[1]
    nt = length(t)
    f_reshape = reshape(f, ny, nt) |> transpose
    return t, f_reshape
end

function get_sensitivity(sol::Solution)
    @unpack t, S, dimensions, coefficients = sol
    if isempty(S)
        error("Sensitivity coefficients S = $S are empty, set save_solution = true \
               and sensitivity != NoSensitivity()")
    end
    ny = dimensions[1]
    np = coefficients[1]
    nt = length(t)
    S_reshape = reshape(S, ny*np, nt) |> transpose
    return t, S_reshape
end

function get_eigenmax(sol::Solution)
    @unpack t, lambda_LR = sol
    if isempty(lambda_LR)
        # TODO: say what options to change
        error("Max eigenvalues lambda_LR = $lambda_LR are empty")
    end
    return t, lambda_LR
end

function get_dimensions(sol::Solution)
    @unpack t, dimensions, coefficients = sol
    return (; nt = length(t), ny = dimensions[1], np = coefficients[1])
end