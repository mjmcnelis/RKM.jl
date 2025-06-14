
function compute_stats!(sol::Solution, save_solution::Bool, adaptive::AdaptiveTimeStep,
                        interpolator::Interpolator, timer::TimeLimit,
                        stage_finder::StageFinder, sensitivity::SensitivityMethod,
                        loop_stats::NamedTuple, config_bytes::Int64)

    @unpack state_jacobian = stage_finder
    JE_y = state_jacobian.evaluations[1]
    # TODO: simplify this
    JE_p = 0
    if !(sensitivity isa NoSensitivity)
        @unpack param_jacobian = sensitivity
        JE_p = param_jacobian.evaluations[1]
    end

    @unpack t, y, f, dy, S, time_steps_taken, JE, rejection_rate, runtime,
            solution_size, sensitivity_size, config_memory, excess_memory = sol

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
                                 stage_finder, save_time)
    @unpack f, y, J, res = update_cache
    @unpack root_method, state_jacobian = stage_finder

    @unpack nt, ny, np = get_dimensions(sol)

    n_samples = round(Int64, nt/10)
    t_idxs = round.(Int64, LinRange(2, nt, n_samples))

    FE_time = 0.0           # functional evaluation time
    JE_time = 0.0           # jacobian evaluation time
    LS_time = 0.0           # linear solve time

    for n in t_idxs
        t = sol.t[n]
        y .= view(sol.y, 1+(n-1)*ny:n*ny)

        FE_stat = @timed ode_wrap!(f, t, y)
        FE_time += FE_stat.time

        if !isempty(J)
            set_wrapper!(ode_wrap!, t)
            JE_stat = @timed evaluate_jacobian!(state_jacobian, J, ode_wrap!, y, f)
            JE_time += JE_stat.time
        end

        # note: linear solve estimate assumes Backward Euler
        if !isempty(J) && !isempty(res) && root_method isa Newton
            y_prev = view(sol.y, 1+(n-2)*ny:(n-1)*ny)
            dt = sol.t[n] - sol.t[n-1]
            @.. res = y - y_prev - dt*f

            LS_stat = @timed begin
                if J isa SparseMatrixCSC                # J <- I - dt.J
                    @.. J.nzval *= (-dt)
                else
                    @.. J *= (-dt)
                end
                for k in diagind(J)
                    J[k] += 1.0
                end
                linear_cache.A = J
                linear_cache.b = res
                solve!(linear_cache)
            end
            LS_time += LS_stat.time
        end
    end

    FE_time *= sol.FE[1] / length(t_idxs)   # TODO: subtract FEs from jacobian
    JE_time *= sol.JE[1] / length(t_idxs)
    LS_time *= sol.JE[1] / length(t_idxs)

    println("")
    println("  Subroutine times (seconds)  ")
    println("---------------------------------")
    println("function evaluations | $(round(FE_time, sigdigits = 4))")
    println("jacobian evaluations | $(round(JE_time, sigdigits = 4))")
    println("linear solve         | $(round(LS_time, sigdigits = 4))")
    println("save solution        | $(round(save_time, sigdigits = 4))")
    println("")

    return nothing
end
