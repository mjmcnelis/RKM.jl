
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
    runtime = sol.runtime
    solution_size = sol.solution_size
    sensitivity_size = sol.sensitivity_size
    config_memory = sol.config_memory
    excess_memory = sol.excess_memory

    time_steps_taken .= timer.total_steps
    # TODO: missing param-jacobian evaluations
    JE .= JE_y + JE_p
    rejection_rate .= compute_step_rejection_rate(adaptive, timer)
    runtime .= loop_stats.time
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
    println("time steps taken     = $(sol.time_steps_taken[1])")
    println("time points saved    = $(length(sol.t))")
    println("step rejection rate  = $(round(sol.rejection_rate[1], sigdigits = 4)) %")
    println("function evaluations = $(sol.FE[1])")
    println("jacobian evaluations = $(sol.JE[1])")
    println("evolution runtime    = $(round(sol.runtime[1], sigdigits = 4)) seconds")
    println("solution size        = $(format_bytes(sol.solution_size[1]))")
    println("sensitivity size     = $(format_bytes(sol.sensitivity_size[1]))")
    println("configuration memory = $(format_bytes(sol.config_memory[1]))")
    println("excess memory        = $(format_bytes(sol.excess_memory[1]))")
end

function get_subroutine_runtimes(ode_wrap!, state_jacobian, root_finder, save_time)

    FE_time = ode_wrap!.subroutine_time
    JE_time = state_jacobian.subroutine_time

    if root_finder isa Newton
        LS_time = root_finder.subroutine_time
    else
        LS_time = [0.0]
    end

    println("")
    println("  Subroutine times (seconds)  ")
    println("---------------------------------")
    println("function evaluations | $(round(FE_time[1], sigdigits = 4))")
    println("jacobian evaluations | $(round(JE_time[1], sigdigits = 4))")
    println("linear solve         | $(round(LS_time[1], sigdigits = 4))")
    println("save solution        | $(round(save_time[1], sigdigits = 4))")
    println("")

    # total_time = round(FE_time[1] + JE_time[1] + LS_time[1] + save_time[1], sigdigits = 4)
    # @show total_time

    return nothing
end
