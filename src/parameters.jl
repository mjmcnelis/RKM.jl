
@kwdef struct Parameters
    adaptive::AdaptiveStepSize
    method::ODEMethod
    t_span::TimeSpan
    timer::TimeLimit
end
