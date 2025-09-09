
function compute_stats!(sol::Solution, save_solution::Bool, adaptive::AdaptiveTimeStep,
                        interpolator::Interpolator, timer::TimeLimit,
                        state_jacobian::JacobianMethod, sensitivity::SensitivityMethod,
                        loop_stats::NamedTuple, config_bytes::Int64)

    JE_y = state_jacobian.evaluations[1]
    # TODO: simplify this
    JE_p = 0
    if !(sensitivity isa NoSensitivity)
        param_jacobian = sensitivity.param_jacobian
        JE_p = param_jacobian.evaluations[1]
    end

    t = sol.t
    y = sol.y
    f = sol.f
    dy = sol.dy
    S = sol.S
    time_steps_taken = sol.time_steps_taken
    JE = sol.JE
    rejection_rate = sol.rejection_rate
    solution_size = sol.solution_size
    sensitivity_size = sol.sensitivity_size
    config_memory = sol.config_memory
    excess_memory = sol.excess_memory

    time_steps_taken .= timer.total_steps
    # TODO: missing param-jacobian evaluations
    JE .= JE_y + JE_p
    rejection_rate .= compute_step_rejection_rate(adaptive, timer)
    solution_size .= sizeof(t) + sizeof(y) + sizeof(f) + sizeof(dy)
    sensitivity_size .= sizeof(S)
    config_memory .= config_bytes
    excess_memory .= loop_stats.bytes
    if !(save_solution && adaptive isa Fixed)
        excess_memory .-= (solution_size .+ sensitivity_size)
    end
    return nothing
end

function get_stats(sol::Solution)
    runtimes = sol.runtimes
    evolution_time = runtimes.evolution_time

    println("time steps taken     = $(sol.time_steps_taken[1])")
    println("time points saved    = $(length(sol.t))")
    println("step rejection rate  = $(round(sol.rejection_rate[1], sigdigits = 4)) %")
    println("function evaluations = $(sol.FE[1])")
    println("jacobian evaluations = $(sol.JE[1])")
    println("evolution runtime    = $(round(evolution_time[1], sigdigits = 4)) seconds")
    println("solution size        = $(format_bytes(sol.solution_size[1]))")
    println("sensitivity size     = $(format_bytes(sol.sensitivity_size[1]))")
    println("configuration memory = $(format_bytes(sol.config_memory[1]))")
    println("excess memory        = $(format_bytes(sol.excess_memory[1]))")
end
