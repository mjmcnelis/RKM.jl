# TODO update docstring
"""
    evolve_ode!(sol::Solution, y0::Union{T, Vector{T}}, t0::T1, tf::Float64,
                dt0::Float64, dy_dt!::Function, options::SolverOptions;
                p::Vector{Float64} = Float64[],
                abstract_params = nothing) where {T <: AbstractFloat,
                                                  T1 <: AbstractFloat}

Required parameters: `sol`, `y0`, `t0`, tf`, `dt0`, `dy_dt!`, `options`
"""
function evolve_ode!(sol::Solution, y0::Union{T, Vector{T}}, t0::T1, tf::Float64,
                     dt0::Float64, dy_dt!::Function, options::SolverOptions,
                     p::Vector{Float64} = Float64[];
                     abstract_params = nothing) where {T <: AbstractFloat,
                                                       T1 <: AbstractFloat}
    config_bytes = @allocated begin
        clear_solution!(sol)
        @unpack precision = options

        # @code_warntype lookup_options(options)
        # q()

        # get solver options
        @unpack adaptive, controller, method, timer, stage_finder,
                interpolator, sensitivity_method, show_progress,
                save_solution, benchmark_subroutines = options

        reset_timer!(timer)

        dimensions = size(y0, 1)
        coefficients = length(p)
        sol.dimensions .= dimensions
        sol.coefficients .= coefficients

        # initial conditions
        y = y0 .|> precision
        y isa Vector ? nothing : y = [y]

        t0 = rationalize(t0) |> precision
        tf = rationalize(tf) |> precision
        dt0 = rationalize(dt0) |> precision

        t  = [t0]
        dt = [dt0, dt0]

        # reconstruction
        method = reconstruct_method(method, precision)
        if !(adaptive isa Fixed)
            adaptive = reconstruct_adaptive(adaptive, method, precision)
            controller = reconstruct_controller(controller, method, adaptive, precision)
        end

        # configure cache
        update_cache = UpdateCache(precision, y, method, adaptive,
                                   dimensions, coefficients,
                                   sensitivity_method, stage_finder)

        @unpack y, y_tmp, f, f_tmp, J, error, S, S_tmp = update_cache

        # create ODE wrapper function
        ode_wrap! = ODEWrapperState([t0], p, abstract_params, dy_dt!)
        # re-using y_tmp for now
        ode_wrap_p! = ODEWrapperParam([t0], y_tmp, abstract_params, dy_dt!) # not used yet

        @unpack linear_method = stage_finder
        # configure linear cache (see src/common.jl in LinearSolve.jl)
        linear_cache = init(LinearProblem(J, error), linear_method;
                            alias_A = true, alias_b = true)

        if method.iteration isa Implicit || !(sensitivity_method isa NoSensitivity)
            stage_finder = set_jacobian_cache(stage_finder, ode_wrap!, f_tmp, y)
        end

        sensitivity_method = set_jacobian_cache(sensitivity_method, ode_wrap_p!,
                                                f_tmp, y, p)

        # for progress meter
        checkpoints = collect(LinRange(t0, tf, 101))[2:end]
        progress = create_progress(100; showspeed = true, color = :gray)
    end

    # sizehint solution
    if save_solution
        sizehint_solution!(adaptive, interpolator, sol, t0, tf, dt0, sensitivity_method)
    end

    # runtime to store solution (S2)
    S2_runtime = 0.0

    # evaluate ODE at initial time and store in f/f_tmp
    ode_wrap!(f, t[1], y)
    @.. f_tmp = f

    loop_stats = @timed begin
        # save initial condition
        if save_solution
            append!(sol.t, t[1])
            append!(sol.y, y)
            if interpolator isa DenseInterpolator
                append!(sol.f, f)
            end
            if !(sensitivity_method isa NoSensitivity)
                append!(sol.S, S)
            end
        end
        # time evolution loop
        while true
            show_progress && monitor_progress(t, progress, checkpoints)
            continue_solver(t, tf, timer) || break

            adjust_final_time_steps!(adaptive, t, dt, tf)

            evolve_one_time_step!(method, adaptive, controller,
                                  t, dt, ode_wrap!, update_cache,
                                  linear_cache, stage_finder,
                                  sensitivity_method, ode_wrap_p!,
                                  interpolator)
            t[1] += dt[1]
            timer.total_steps[1] += 1

            if save_solution
                # TODO: this is pretty cumbersome
                if benchmark_subroutines && length(sol.t) % 10 == 0
                    S2_stat = @timed begin
                        append!(sol.t, t[1])
                        append!(sol.y, y_tmp)
                        if interpolator isa DenseInterpolator
                            append!(sol.f, f_tmp)
                        end
                        if !(sensitivity_method isa NoSensitivity)
                            append!(sol.S, S_tmp)
                        end
                    end
                    S2_runtime += 10.0*S2_stat.time
                else
                    append!(sol.t, t[1])
                    append!(sol.y, y_tmp)
                    if interpolator isa DenseInterpolator
                        append!(sol.f, f_tmp)
                    end
                    if !(sensitivity_method isa NoSensitivity)
                        append!(sol.S, S_tmp)
                    end
                end
            end
            # store updated values in tmp caches
            @.. y = y_tmp
            @.. S = S_tmp
            @.. f = f_tmp
        end
    end

    # TODO: move to compute_stats!
    sol.FE[1] = ode_wrap!.FE[1] + ode_wrap_p!.FE[1]

    compute_stats!(sol, save_solution, adaptive, interpolator, timer,
                   stage_finder, loop_stats, config_bytes)

    if benchmark_subroutines && save_solution
        get_subroutine_runtimes(sol, ode_wrap!, update_cache, linear_cache,
                                stage_finder, S2_runtime, n_samples = 100)
    end

    return nothing
end

"""
    evolve_ode(y0::Union{T, Vector{T}}, t0::T1, tf::Float64, dt0::Float64,
               dy_dt!::Function, options::SolverOptions;
               p::Vector{Float64} = Float64[],
               model_parameters = nothing) where {T <: AbstractFloat,
                                                  T1 <: AbstractFloat}

Required parameters: `y0`, `t0`, `tf`, `dt0`, `dy_dt!`, `options`
"""
function evolve_ode(y0::Union{T, Vector{T}}, t0::T1, tf::Float64, dt0::Float64,
                    dy_dt!::Function, options::SolverOptions,
                    p::Vector{Float64} = Float64[];
                    abstract_params = nothing) where {T <: AbstractFloat,
                                                      T1 <: AbstractFloat}
    sol = Solution(options)
    evolve_ode!(sol, y0, t0, tf, dt0, dy_dt!, options, p; abstract_params)
    return sol
end
