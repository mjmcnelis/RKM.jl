
function compute_stats!(sol::Solution, save_solution::Bool, adaptive::AdaptiveStepSize,
                        interpolator::Interpolator, timer::TimeLimit,
                        stage_finder::StageFinder, loop_stats::NamedTuple,
                        config_bytes::Int64)

    @unpack jacobian_method = stage_finder
    @unpack time_steps_taken, JE, rejection_rate, runtime,
            solution_size, config_memory, excess_memory = sol

    time_steps_taken .= timer.total_steps
    JE .= jacobian_method.evaluations
    rejection_rate .= compute_step_rejection_rate(adaptive, timer)
    runtime .= loop_stats.time
    solution_size .= sizeof(sol.y) + sizeof(sol.t)
    config_memory .= config_bytes
    excess_memory .= loop_stats.bytes
    # I sizehint if Fixed (any interpolator) or use
    if !(save_solution && (adaptive isa Fixed || interpolator isa DenseInterpolator))
        excess_memory .-= solution_size
    end
    return nothing
end

function get_stats(sol::Solution)
    @unpack t, time_steps_taken, rejection_rate, FE, JE, runtime,
            solution_size, config_memory, excess_memory = sol
    println("time steps taken     = $(time_steps_taken[1])")
    println("time points saved    = $(length(t))")
    println("step rejection rate  = $(round(rejection_rate[1], sigdigits = 4)) %")
    println("function evaluations = $(FE[1])")
    println("jacobian evaluations = $(JE[1])")
    println("evolution runtime    = $(round(runtime[1], sigdigits = 4)) seconds")
    println("solution size        = $(format_bytes(solution_size[1]))")
    println("config memory        = $(format_bytes(config_memory[1]))")
    println("excess memory        = $(format_bytes(excess_memory[1]))")
end
