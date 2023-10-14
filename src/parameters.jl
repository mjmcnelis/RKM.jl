"""
Parameters for the ODE solver.
"""
@kwdef struct Parameters
    """Adaptive time step method"""
    adaptive::AdaptiveStepSize
    """ODE solver method"""
    method::ODEMethod
    """Time range of ODE evolution"""
    t_range::TimeRange
    """Timer for ODE solver"""
    timer::TimeLimit = TimeLimit(; wtime_min = 60)
    """Adaptive time step controller"""
    controller::Controller = PIDControlK()
    """Stage finder for implicit ODE methods"""
    stage_finder::StageFinder = ImplicitStageFinder()
end
# note: need to use abstract types (e.g. ::StageFinder)
#       to avoid excess allocations in evolve loop

function Parameters(dict::Dict)
    return Parameters(; (Symbol(k) => v for (k,v) in dict)...)
end