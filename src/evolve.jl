# TODO update docstring
"""
    evolve_ode!(sol::Solution{T1}, y0::Vector{T}, t0::T, tf::Float64,
                dt0::Float64, dy_dt!::Function, options::SolverOptions{T1},
                p::Vector{Float64} = Float64[];
                abstract_params = nothing) where {T <: AbstractFloat,
                                                  T1 <: AbstractFloat}

Required parameters: `sol`, `y0`, `t0`, tf`, `dt0`, `dy_dt!`, `options`
"""
function evolve_ode!(sol::Solution{T1}, y0::Vector{T}, t0::T, tf::Float64,
                     dt0::Float64, dy_dt!::Function, options::SolverOptions{T1},
                     p::Vector{Float64} = Float64[];
                     abstract_params = nothing) where {T <: AbstractFloat,
                                                       T1 <: AbstractFloat}
    config_bytes = @allocated begin
        clear_solution!(sol)

        # get solver options
        @unpack adaptive, method, timer, stage_finder, interpolator,
                sensitivity, show_progress, save_solution,
                save_time_derivative, benchmark_subroutines, precision = options

        reset_timer!(timer)

        dimensions = size(y0, 1)
        coefficients = length(p)
        sol.dimensions .= dimensions
        sol.coefficients .= coefficients

        # initial conditions
        y = y0 .|> precision

        t0 = rationalize(t0) |> precision
        tf = rationalize(tf) |> precision
        dt0 = rationalize(dt0) |> precision

        t  = [t0]
        dt = [dt0, dt0]

        # reconstruct: method, adaptive, stage_finder, sensitivity
        # can these be done in SolverOptions constructor?

        # reconstruction
        method = reconstruct_method(method, precision)

        # TODO: put this before method reconstruction?
        #       or maybe unpack and pass local order
        adaptive = reconstruct_adaptive(adaptive, method)

        # configure cache
        update_cache = UpdateCache(precision, y, method, adaptive, dimensions,
                                   coefficients, sensitivity, stage_finder)

        @unpack y, y_tmp, f, f_tmp, dy, J, res, S, S_tmp, lambda_LR, x0 = update_cache

        # create ODE wrappers
        ode_wrap_y! = ODEWrapperState([t0], p, abstract_params, dy_dt!)#, method)
        ode_wrap_p! = ODEWrapperParam([t0], y_tmp, abstract_params, dy_dt!)#, method)

        @unpack iteration = method
        if iteration isa Implicit || !(sensitivity isa NoSensitivity)
            stage_finder = reconstruct_stage_finder(stage_finder, ode_wrap_y!, f_tmp, y)
        end

        @unpack linear_method, eigenmax, state_jacobian = stage_finder
        # configure linear cache (see src/common.jl in LinearSolve.jl)
        linear_cache = init(LinearProblem(J, res), linear_method;
                            alias_A = true, alias_b = true)

        sensitivity = reconstruct_sensitivity(sensitivity, ode_wrap_p!, f_tmp, p)

        # for progress meter
        checkpoints = collect(LinRange(t0, tf, 101))[2:end]
        progress = Progress(100; desc = "  Progress: ", showspeed = true, color = :gray)

        # time to store solution
        save_time = [0.0]

        # evaluate ODE at initial time and store in f/f_tmp
        ode_wrap_y!(f, t[1], y)
        @.. f_tmp = f

        # initialize eigenvalue/eigenvector
        @.. y_tmp = y
        compute_max_eigenvalue!(eigenmax, lambda_LR, x0, J, state_jacobian,
                                ode_wrap_y!, y_tmp, f_tmp)
    end

    # sizehint solution
    if save_solution
        @unpack stages = method
        sizehint_solution!(adaptive, interpolator, sol, t0, tf, dt0, sensitivity,
                           save_time_derivative, stages, iteration, eigenmax)
    end

    loop_stats = @timed begin
        # save initial condition
        if save_solution
            # TODO: this output function still allocates
            # initial_output!(sol, update_cache, t, options)
            append!(sol.t, t[1])
            append!(sol.y, y)
            if save_time_derivative || interpolator isa CubicHermite
                append!(sol.f, f)
            end
            if !(sensitivity isa NoSensitivity)
                append!(sol.S, S)
            end
            if iteration isa Implicit && !(eigenmax isa NoEigenMax)
                append!(sol.lambda_LR, lambda_LR)
            end
        end
        # time evolution loop
        while true
            if show_progress
                set_current_system_time!(timer)
                monitor_progress(progress, checkpoints, t, timer, dt)
            end
            if !continue_solver(t, dt, tf, timer, show_progress)
                break
            end
            adjust_final_time_steps!(t, dt, tf)

            evolve_one_time_step!(method, adaptive, t, dt, ode_wrap_y!,
                                  update_cache, linear_cache, stage_finder,
                                  sensitivity, ode_wrap_p!, interpolator)
            t[1] += dt[1]
            timer.total_steps[1] += 1

            # note: missing dy evaluation at final time (not critical but awkward)
            if save_solution
                output_solution!(sol, save_time, update_cache, t, options, timer)
                # TODO: move to output_solution
                if iteration isa Implicit && !(eigenmax isa NoEigenMax)
                    append!(sol.lambda_LR, lambda_LR)
                end
            end
            # get updated values from tmp caches
            @.. y = y_tmp
            @.. S = S_tmp
            @.. f = f_tmp
        end
    end

    # TODO: move to compute_stats!
    sol.FE[1] = ode_wrap_y!.FE[1] + ode_wrap_p!.FE[1]

    compute_stats!(sol, save_solution, adaptive, interpolator, timer,
                   stage_finder, sensitivity, loop_stats, config_bytes)

    if benchmark_subroutines && save_solution
        get_subroutine_runtimes(sol, ode_wrap_y!, update_cache, linear_cache,
                                stage_finder, save_time[1])
    end

    return nothing
end

"""
    evolve_ode(y0::Vector{T}, t0::T, tf::Float64, dt0::Float64,
               dy_dt!::Function, options::SolverOptions{T1},
               p::Vector{Float64} = Float64[];
               abstract_params = nothing) where {T <: AbstractFloat,
                                                    T1 <: AbstractFloat}

Required parameters: `y0`, `t0`, `tf`, `dt0`, `dy_dt!`, `options`
"""
function evolve_ode(y0::Vector{T}, t0::T, tf::Float64, dt0::Float64,
                    dy_dt!::Function, options::SolverOptions{T1},
                    p::Vector{Float64} = Float64[];
                    abstract_params = nothing) where {T <: AbstractFloat,
                                                      T1 <: AbstractFloat}
    sol = Solution(options)
    evolve_ode!(sol, y0, t0, tf, dt0, dy_dt!, options, p; abstract_params)
    return sol
end
