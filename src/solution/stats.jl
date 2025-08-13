
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
    println("config memory        = $(format_bytes(sol.config_memory[1]))")
    println("excess memory        = $(format_bytes(sol.excess_memory[1]))")
end

function get_subroutine_runtimes(sol, ode_wrap!, update_cache, root_finder,
                                 state_jacobian, save_time)

    f = update_cache.f
    y = update_cache.y
    J = update_cache.J
    res = update_cache.res

    nt = length(sol.t)
    ny = sol.dimensions[1]

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
        if !isempty(J) && !isempty(res) && root_finder isa Newton
            linear_cache = root_finder.linear_cache

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
