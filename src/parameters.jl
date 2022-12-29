"""
Parameters for the ODE solver.
"""
@kwdef struct Parameters
    """Adaptive time step method"""
    adaptive::AdaptiveStepSize
    """ODE solver method"""
    method::ODEMethod
    """Time span of ODE evolution"""
    t_span::TimeSpan
    """Timer for ODE solver"""
    timer::TimeLimit = TimeLimit(; wtime_min = 60, frequency = 100)
end
