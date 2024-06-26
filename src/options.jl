"""
SolverOptions for the ODE solver.
"""
@kwdef struct SolverOptions
    """Adaptive time step method"""
    adaptive::AdaptiveStepSize
    """ODE solver method"""
    method::ODEMethod
    """Timer for ODE solver"""
    timer::TimeLimit = TimeLimit(; wtime_min = Inf)
    """Adaptive time step controller"""
    controller::Controller = TimeStepController()
    """Stage finder for implicit ODE methods"""
    stage_finder::StageFinder = ImplicitStageFinder()
    """Interpolation method for dense output"""
    interpolator::Interpolator = NoInterpolator()
    """Determines whether or not the solution is stored"""
    save_solution::Bool = true
    """Determines whether or not the progress meter is displayed"""
    show_progress::Bool = false
end
# note: need to use abstract types (e.g. ::StageFinder)
#       to avoid excess allocations in evolve loop

function SolverOptions(dict::Dict)
    return SolverOptions(; (Symbol(k) => v for (k,v) in dict)...)
end

function lookup_options(options::SolverOptions)
    @unpack adaptive, method, timer, controller, stage_finder,
            interpolator, save_solution, show_progress = options
    return nothing
end