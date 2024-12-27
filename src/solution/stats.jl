
function compute_stats!(sol::Solution, save_solution::Bool, adaptive::AdaptiveStepSize,
                        interpolator::Interpolator, timer::TimeLimit,
                        stage_finder::StageFinder, loop_stats::NamedTuple,
                        config_bytes::Int64)

    @unpack jacobian_method = stage_finder
    @unpack time_steps_taken, JE, rejection_rate, runtime,
            solution_size, sensitivity_size, config_memory, excess_memory = sol

    time_steps_taken .= timer.total_steps
    JE .= jacobian_method.evaluations
    rejection_rate .= compute_step_rejection_rate(adaptive, timer)
    runtime .= loop_stats.time
    solution_size .= sizeof(sol.y) + sizeof(sol.t)
    sensitivity_size .= sizeof(sol.S)
    config_memory .= config_bytes
    excess_memory .= loop_stats.bytes
    # I sizehint if Fixed (any interpolator) or use
    if !(save_solution && (adaptive isa Fixed || interpolator isa DenseInterpolator))
        excess_memory .-= (solution_size .+ sensitivity_size)
    end
    return nothing
end

function get_stats(sol::Solution)
    @unpack t, time_steps_taken, rejection_rate, FE, JE, runtime,
            solution_size, sensitivity_size, config_memory, excess_memory = sol
    println("time steps taken     = $(time_steps_taken[1])")
    println("time points saved    = $(length(t))")
    println("step rejection rate  = $(round(rejection_rate[1], sigdigits = 4)) %")
    println("function evaluations = $(FE[1])")
    println("jacobian evaluations = $(JE[1])")
    println("evolution runtime    = $(round(runtime[1], sigdigits = 4)) seconds")
    println("solution size        = $(format_bytes(solution_size[1]))")
    println("sensitivity size     = $(format_bytes(sensitivity_size[1]))")
    println("config memory        = $(format_bytes(config_memory[1]))")
    println("excess memory        = $(format_bytes(excess_memory[1]))")
end

function get_subroutine_runtimes(sol, ode_wrap!, update_cache, linear_cache,
                                 stage_finder; n_samples = 100)
    @unpack f, y, J, error = update_cache
    @unpack root_method, jacobian_method = stage_finder

    ny = sol.dimensions[1]
    nt = length(sol.t)
    FE_dummy = [0]          # placeholder variable

    # sample subset of time indices (could try StatsBase.sample)
    t_idxs = round.(Int64, LinRange(2, nt, min(nt, n_samples)))

    # sample subset of state variables
    y_sol = typeof(sol.y)()
    @time sizehint!(y_sol, ny*length(t_idxs))

    FE_runtime = 0.0        # functional evaluation time
    JE_runtime = 0.0        # jacobian evaluation time
    LS_runtime = 0.0        # linear solve time
    S2_runtime = 0.0        # save solution time

    for n in t_idxs
        t = sol.t[n]
        y .= view(sol.y, 1+(n-1)*ny:n*ny)

        S2_stat = @timed append!(y_sol, y)
        S2_runtime += S2_stat.time

        ode_wrap!.t[1] = t

        FE_stat = @timed ode_wrap!(f, t, y)
        FE_runtime += FE_stat.time

        if !isempty(J) && root_method isa Newton
            JE_stat = @timed evaluate_system_jacobian!(jacobian_method, FE_dummy,
                                                    J, ode_wrap!, y, f)
            JE_runtime += JE_stat.time

            # note: linear solve estimate assumes Backward Euler
            LS_stat = @timed begin
                y_prev = view(sol.y, 1+(n-2)*ny:(n-1)*ny)
                dt = sol.t[n] - sol.t[n-1]

                @.. J *= dt
                for k in diagind(J)
                    J[k] += 1.0
                end
                @.. error = y - y_prev - dt*f
                linear_cache.A = J
                linear_cache.b = error
                solve!(linear_cache)
            end
            LS_runtime += LS_stat.time
        end
    end

    FE_runtime *= sol.FE[1] / length(t_idxs)    # TODO: subtract FEs from jacobian
    JE_runtime *= sol.JE[1] / length(t_idxs)
    LS_runtime *= sol.JE[1] / length(t_idxs)
    S2_runtime *= nt / length(t_idxs)

    println("")
    println("  Subroutine runtimes (seconds)  ")
    println("---------------------------------")
    println("function evaluations | $(round(FE_runtime, sigdigits = 4))")
    println("jacobian evaluations | $(round(JE_runtime, sigdigits = 4))")
    println("linear solve         | $(round(LS_runtime, sigdigits = 4))")
    println("save solution        | $(round(S2_runtime, sigdigits = 4))")
    println("")

    return nothing
end
