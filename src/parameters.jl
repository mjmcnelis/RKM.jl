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
    timer::TimeLimit
end
