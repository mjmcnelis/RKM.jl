
function compute_stats!(sol::Solution{T}, save_solution::Bool, adaptive::AdaptiveTimeStep,
                        timer::TimeLimit, config::RKMConfig, loop_stats::NamedTuple,
                        config_bytes::Int64) where T <: AbstractFloat

    ode_wrap_y! = config.ode_wrap_y!
    ode_wrap_p! = config.ode_wrap_p!
    state_jacobian = config.state_jacobian
    sensitivity = config.sensitivity

    t = sol.t
    y = sol.y
    f = sol.f
    dy = sol.dy
    S = sol.S
    total_steps = sol.total_steps
    FE = sol.FE
    JE = sol.JE
    rejection_rate = sol.rejection_rate
    solution_size = sol.solution_size
    sensitivity_size = sol.sensitivity_size
    config_memory = sol.config_memory
    excess_memory = sol.excess_memory

    total_steps[1] = timer.total_steps[1]
    FE[1] = sum(ode_wrap_y!.evaluations) + ode_wrap_p!.evaluations[1]

    JE_y = state_jacobian.evaluations[1]
    if sensitivity isa NoSensitivity
        JE_p = 0
    else
        JE_p = sensitivity.param_jacobian.evaluations[1]
    end
    JE[1] = JE_y[1] + JE_p[1]

    rejection_rate[1] = compute_step_rejection_rate(adaptive, timer)
    solution_size[1] = sum(x -> sizeof(x), (t, y, f, dy))
    sensitivity_size[1] = sizeof(S)

    # correct memory if use BigFloat precision
    if T == BigFloat
        n_sol = sum(x -> length(x), (t, y, f, dy))
        n_sen = length(S)
        BF_bytes = summarysize(BigFloat(1.0))

        solution_size[1] += n_sol*BF_bytes
        sensitivity_size[1] += n_sen*BF_bytes
    end

    # @show summarysize(sol.t) + summarysize(sol.y) solution_size[1]
    # @show summarysize(sol.S) sensitivity_size[1]

    config_memory[1] = config_bytes
    excess_memory[1] = loop_stats.bytes
    if !(save_solution && adaptive isa Fixed)
        excess_memory[1] -= (solution_size[1] + sensitivity_size[1])
    end

    return nothing
end

function get_stats(sol::Solution)
    runtimes = sol.runtimes
    evolution_time = runtimes.evolution_time

    println("total steps          = $(sol.total_steps[1])")
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
