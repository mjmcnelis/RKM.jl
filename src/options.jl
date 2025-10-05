"""
$(TYPEDEF)

List of options for the ODE solver.

$(TYPEDFIELDS)
"""
@kwdef struct SolverOptions{T <: AbstractFloat}
    """ODE solver method"""
    method::ODEMethod
    """Adaptive time step method"""
    adaptive::AdaptiveTimeStep
    """Timer for ODE solver"""
    timer::TimeLimit = TimeLimit()
    """State Jacobian method for implicit ODE solvers or sensitivity analysis"""
    state_jacobian::JacobianMethod = FiniteJacobian()
    """Root finder for implicit ODE solvers"""
    root_finder::RootFinderMethod = Newton()
    """Max eigenvalue method for implicit ODE solvers"""
    eigenmax::EigenMaxMethod = NoEigenMax()
    """Sensitivity method"""
    sensitivity::SensitivityMethod = NoSensitivity()
    """Interpolation method for dense output"""
    interpolator::Interpolator = NoInterpolation()
    """Store numerical solution (time and state variables)"""
    save_solution::Bool = true
    """Additionally output time derivatives"""
    save_time_derivative::Bool = false
    """Display progress meter"""
    show_progress::Bool = false
    """Estimate runtime of core subroutines (e.g. function evaluations)"""
    time_subroutine::Bool = false
    """Numerical precision of the solver and solution"""
    precision::Type{T} = Float64
end
# note: need to use abstract types (e.g. ::StageFinder)
#       to avoid excess allocations in evolve loop

function SolverOptions(dict::Dict)
    return SolverOptions(; (Symbol(k) => v for (k,v) in dict)...)
end

function lookup_options(options::SolverOptions)
    adaptive = options.adaptive
    method = options.method
    timer = options.timer
    state_jacobian = options.state_jacobian
    root_finder = options.root_finder
    eigenmax = options.eigenmax
    sensitivity = options.sensitivity
    interpolator = options.interpolator
    save_solution = options.save_solution
    save_time_derivative = options.save_time_derivative
    show_progress = options.show_progress
    time_subroutine = options.time_subroutine
    precision = options.precision
    return nothing
end

"""
    Solution(options::SolverOptions)

Outer constructor for `Solution`.

Required parameters: `options`
"""
function Solution(options::SolverOptions)
    return Solution(options.precision)
end