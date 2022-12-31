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
    timer::TimeLimit = TimeLimit(; wtime_min = 60, frequency = 100)
end
