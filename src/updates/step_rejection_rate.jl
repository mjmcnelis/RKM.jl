
_step_rejection_rate(::Iteration, ::AdaptiveStepSize, args...) = 0.0

function _step_rejection_rate(::Explicit, ::Doubling, stages, steps, FE, args...)
    # reconstruct total number of attempts in step doubling routine
    attempts = (FE - steps) / (3*stages - 2)
    return 100.0 * (1.0 - steps/attempts)
end

function _step_rejection_rate(::Explicit, ::Embedded, stages, steps, FE, fsal)
    fsal_stage = fsal isa FSAL ? 1 : 0
    # reconstruct total number of attempts in embedded routine
    attempts = (FE - steps*(1 - fsal_stage)) / (stages - 1)
    return 100.0 * (1.0 - steps/attempts)
end

function compute_step_rejection_rate!(sol::Solution, method::ODEMethod,
                                      adaptive::AdaptiveStepSize, timer::TimeLimit)
    steps = timer.counter[1] - 1
    FE = sol.FE[1]
    @unpack stages, fsal, iteration = method
    @unpack rejection_rate = sol

    rejection_rate .= _step_rejection_rate(iteration, adaptive, stages, steps, FE, fsal)
    return nothing
end