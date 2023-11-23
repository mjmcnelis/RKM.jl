
function compute_stats!(sol::Solution, options, stage_finder::StageFinder,
                        loop_stats::NamedTuple, config_bytes::Int64)

    @unpack save_solution, adaptive = options
    @unpack jacobian_method = stage_finder
    @unpack JE, runtime, solution_size, config_memory, excess_memory = sol

    # TODO: fix step rejection rate
    # compute_step_rejection_rate!(sol, method, adaptive, timer)
    JE .= jacobian_method.evaluations
    runtime .= loop_stats.time
    solution_size .= sizeof(sol.y) + sizeof(sol.t)
    config_memory .= config_bytes
    excess_memory .= loop_stats.bytes
    if !(save_solution && adaptive isa Fixed)
        excess_memory .-= solution_size
    end
    return  nothing
end

function get_stats(sol::Solution)
    @unpack y, t, FE, JE, rejection_rate, runtime,
            solution_size, config_memory, excess_memory = sol
    println("time steps           = $(length(t))")
    println("step rejection rate  = $(sol.rejection_rate[1]) %")
    println("function evaluations = $(FE[1])")
    println("jacobian evaluations = $(JE[1])")
    println("evolution runtime    = $(round(runtime[1], sigdigits = 4)) seconds")
    println("solution size        = $(format_bytes(solution_size[1]))")
    println("config memory        = $(format_bytes(config_memory[1]))")
    println("excess memory        = $(format_bytes(excess_memory[1]))")
end
