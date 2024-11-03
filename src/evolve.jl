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
        @unpack FE#=, JE=# = sol
        @unpack precision = options

        # @code_warntype lookup_options(options)
        # q()

        # get solver options
        @unpack adaptive, controller, method, timer, stage_finder,
                interpolator, sensitivity_method, show_progress,
                save_solution = options

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

        t  = [t0, t0]
        dt = [dt0, dt0]

        # reconstruction
        method = reconstruct_method(method, precision)
        interpolator = reconstruct_interpolator(interpolator, t0, tf)
        if !(adaptive isa Fixed)
            adaptive = reconstruct_adaptive(adaptive, method, precision)
            controller = reconstruct_controller(controller, method, adaptive, precision)
        end

        # configure cache
        update_cache = UpdateCache(precision, y, method, adaptive,
                                   dimensions, coefficients,
                                   sensitivity_method)

        @unpack y, y_tmp, f, f_tmp, J, error, S, S_tmp = update_cache

        # create ODE wrapper function
        ode_wrap! = ODEWrapperState([t0], p, abstract_params, dy_dt!)
        # re-using y_tmp for now
        ode_wrap_p! = ODEWrapperParam([t0], y_tmp, abstract_params, dy_dt!) # not used yet

        J[diagind(J)] .= 1.0
        # TODO: have option to use sparse jacobian
        # J = sparse(J)

        @unpack linear_method = stage_finder
        # configure linear cache (see src/common.jl in LinearSolve.jl)
        # note: assumes LU factorization for now
        linear_cache = init(LinearProblem(J, error), linear_method;
                            alias_A = true, alias_b = true)

        if method.iteration isa Implicit || !(sensitivity_method isa NoSensitivity)
            stage_finder = set_jacobian_cache(stage_finder, ode_wrap!, f_tmp, y)
        end

        sensitivity_method = set_jacobian_cache(sensitivity_method, p, y)

        # for progress meter
        checkpoints = collect(LinRange(t0, tf, 101))[2:end]
        progress = create_progress(100; showspeed = true, color = :gray)
    end

    # sizehint solution
    if save_solution
        sizehint_solution!(adaptive, interpolator, sol, t0, tf, dt0, sensitivity_method)
    end

    loop_stats = @timed begin
        # save initial condition
        if save_solution
            append!(sol.y, y)
            append!(sol.S, S)
            append!(sol.t, t[1])
        end
        # time evolution loop
        while true
            show_progress && monitor_progress(t, progress, checkpoints)
            continue_solver(t, tf, timer) || break

            evolve_one_time_step!(method, adaptive, controller,
                                  FE, t, dt, ode_wrap!, update_cache,
                                  linear_cache, stage_finder,
                                  sensitivity_method, ode_wrap_p!)
            t[2] = t[1]
            t[1] += dt[1]
            timer.total_steps[1] += 1

            if save_solution
                interpolate_solution!(interpolator, sol, update_cache, t)
            end
            # note: make sure updated value is stored in y_tmp, S_tmp
            @.. y = y_tmp
            @.. S = S_tmp
        end
    end
    compute_stats!(sol, save_solution, adaptive, interpolator, timer,
                   stage_finder, loop_stats, config_bytes)
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
