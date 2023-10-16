"""
Parameters for the ODE solver.
"""
@kwdef struct Parameters
    """Adaptive time step method"""
    adaptive::AdaptiveStepSize
    """ODE solver method"""
    method::ODEMethod
    """Time range of ODE evolution"""
    t_range::TimeRange = TimeRange(; t0 = 0.0, tf = 1.0)
    """Timer for ODE solver"""
    timer::TimeLimit = TimeLimit(; wtime_min = 60)
    """Adaptive time step controller"""
    controller::Controller = TimeStepController()
    """Stage finder for implicit ODE methods"""
    stage_finder::StageFinder = ImplicitStageFinder()
    """Determines whether or not the solution is stored"""
    save_solution::Bool = true
    """Determines whether or not the progress meter is displayed"""
    show_progress::Bool = false
    """Determines whether or not static array types are used"""
    static_array::Bool = false
end
# note: need to use abstract types (e.g. ::StageFinder)
#       to avoid excess allocations in evolve loop

function Parameters(dict::Dict)
    return Parameters(; (Symbol(k) => v for (k,v) in dict)...)
end