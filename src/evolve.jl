# TODO update docstring
"""
    evolve_ode!(sol::Solution, y0::Union{T, Vector{T}}, t0::T1, tf::Float64,
                dt0::Float64, dy_dt!::Function, options::SolverOptions;
                model_parameters = nothing) where {T <: AbstractFloat,
                                                   T1 <: AbstractFloat}

Required parameters: `sol`, `y0`, `t0`, tf`, `dt0`, `dy_dt!`, `options`
"""
function evolve_ode!(sol::Solution, y0::Union{T, Vector{T}}, t0::T1, tf::Float64,
                     dt0::Float64, dy_dt!::Function, options::SolverOptions;
                     model_parameters = nothing) where {T <: AbstractFloat,
                                                        T1 <: AbstractFloat}
    config_bytes = @allocated begin
        clear_solution!(sol)
        @unpack precision, FE#=, JE=# = sol

        # @code_warntype lookup_options(options)
        # q()

        @unpack adaptive, controller, method, timer, stage_finder,
                interpolator, show_progress, save_solution = options

        reset_timer!(timer)
        reset_attempts!(adaptive)

        method = reconstruct_method(method, precision)
        if !(adaptive isa Fixed)
            controller = reconstruct_controller(controller, method, adaptive, precision)
        end

        dimensions = size(y0, 1)
        sol.dimensions .= dimensions

        # initial conditions
        y = y0 .|> precision
        y isa Vector ? nothing : y = [y]

        t0 = rationalize(t0) |> precision
        tf = rationalize(tf) |> precision
        dt0 = rationalize(dt0) |> precision

        t  = [t0, t0]
        dt = [dt0, dt0]

        # create ODE wrapper function
        ode_wrap! = ODEWrapper([t0], model_parameters, dy_dt!)

        # configure cache
        update_cache = UpdateCache(precision, y, method, adaptive, dimensions)

        @unpack y, f, f_tmp, J, error = update_cache

        J[diagind(J)] .= 1.0
        # TODO: have option to use sparse jacobian
        # J = sparse(J)
        fact = lu(J)

        # configure linear cache (see src/common.jl in LinearSolve.jl)
        # note: assumes LU factorization for now
        linear_cache = init(LinearProblem(J, error), LUFactorization())#; alias_A = false, alias_b = false)
        linear_cache.cacheval = fact
        linear_cache.isfresh = false

        if method.iteration isa Implicit
            stage_finder = set_jacobian_cache(stage_finder, ode_wrap!, f_tmp, y)
        end

        # for progress meter
        checkpoints = collect(LinRange(t0, tf, 101))[2:end]
        progress = Progress(100; showspeed = true, color = :gray)
    end

    if save_solution && adaptive isa Fixed
        sizehint_solution!(sol, t0, tf, dt0, dimensions)
    end

    # TODO: look into @code_warntype
    loop_stats = @timed begin
        # save initial condition
        if save_solution
            append!(sol.y, y)
            append!(sol.t, t[1])
        end
        # time evolution loop
        while true
            show_progress && monitor_progess(t, progress, checkpoints)
            continue_solver(t, tf, timer) || break

            evolve_one_time_step!(method, adaptive, controller,
                                  FE, t, dt, ode_wrap!, update_cache,
                                  linear_cache, stage_finder)
            t[2] = t[1]
            t[1] += dt[1]
            timer.total_steps[1] += 1

            if save_solution
                interpolate_solution!(interpolator, sol, update_cache, t, t0, tf)
            end
        end
    end
    compute_stats!(sol, save_solution, adaptive, timer,
                   stage_finder, loop_stats, config_bytes)
    return nothing
end

"""
    evolve_ode(y0::Union{T, Vector{T}}, t0::T1, tf::Float64, dt0::Float64,
               dy_dt!::Function, options::SolverOptions,
               precision::Type{T2} = Float64;
               model_parameters = nothing) where {T <: AbstractFloat,
                                                  T1 <: AbstractFloat,
                                                  T2 <: AbstractFloat}

Required parameters: `y0`, `t0`, `tf`, `dt0`, `dy_dt!`, `options`, `precision`
"""
function evolve_ode(y0::Union{T, Vector{T}}, t0::T1, tf::Float64, dt0::Float64,
                    dy_dt!::Function, options::SolverOptions,
                    precision::Type{T2} = Float64;
                    model_parameters = nothing) where {T <: AbstractFloat,
                                                       T1 <: AbstractFloat,
                                                       T2 <: AbstractFloat}
    sol = Solution(precision)
    evolve_ode!(sol, y0, t0, tf, dt0, dy_dt!, options; model_parameters)
    return sol
end
