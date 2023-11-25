# TODO update docstring
"""
    evolve_ode!(sol::Solution, y0::Union{T, Vector{T}}, t0::T1, tf::Float64,
                dt0::Float64, dy_dt!::Function, options::SolverOptions;
                model_parameters = nothing) where {T <: AbstractFloat,
                                                   T1 <: AbstractFloat}

Required parameters: `sol`, `y0`, `dt0`, `dy_dt!`, `options`
"""
function evolve_ode!(sol::Solution, y0::Union{T, Vector{T}}, t0::T1, tf::Float64,
                     dt0::Float64, dy_dt!::Function, options::SolverOptions;
                     model_parameters = nothing) where {T <: AbstractFloat,
                                                        T1 <: AbstractFloat}
    config_bytes = @allocated begin
        clear_solution!(sol)
        @unpack precision, FE#=, JE=# = sol

        @unpack adaptive, controller, method, t_range, timer,
                stage_finder, static_array, show_progress, save_solution = options

        @set! t_range.t0 = t0
        @set! t_range.tf = tf
        @unpack t0, tf = t_range

        reset_timer!(timer)
        reset_attempts!(adaptive)

        method = reconstruct_method(method, precision)
        if !(adaptive isa Fixed)
            controller = reconstruct_controller(controller, method, adaptive, precision)
        end

        @unpack stages = method

        dimensions = size(y0, 1)
        sol.dimensions .= dimensions

        # initial conditions
        y = y0 .|> precision
        y isa Vector ? nothing : y = [y]
        static_array ? y = MVector{dimensions}(y...) : nothing

        t0 = rationalize(t0) |> precision
        tf = rationalize(tf) |> precision
        dt0 = rationalize(dt0) |> precision

        t  = precision == BigFloat ? [t0] : MVector{1}(t0)
        dt = precision == BigFloat ? [dt0, dt0] : MVector{2}(dt0, dt0)

        # create ODE wrapper function
        ode_wrap! = ODEWrapper(deepcopy(t), model_parameters, dy_dt!)

        update_cache = if static_array
            StaticUpdateCache(; method, adaptive, precision, dimensions, stages)
        else
            UpdateCache(; method, adaptive, precision, dimensions, stages)
        end

        @unpack J, f, f_tmp = update_cache

        J[diagind(J)] .= 1.0
        # TODO: have option to use sparse jacobian
        # J = sparse(J)
        fact = lu(J)

        # configure linear cache (see src/common.jl in LinearSolve.jl)
        # note: assumes LU factorization for now
        linear_cache = init(LinearProblem(J, f), LUFactorization())#; alias_A = false, alias_b = false)
        linear_cache.cacheval = fact
        linear_cache.isfresh = false

        stage_finder = set_jacobian_cache(stage_finder, ode_wrap!, f_tmp, y)

        # for progress meter
        checkpoints = collect(LinRange(t0, tf, 101))[2:end]
        progress = Progress(100; showspeed = true, color = :gray)
    end

    if save_solution && adaptive isa Fixed
        sizehint_solution!(sol, t_range, dt0, dimensions)
    end

    # TODO: look into @code_warntype
    loop_stats = @timed while true
        if save_solution
            append!(sol.y, y)
            append!(sol.t, t[1])
        end

        show_progress && monitor_progess(t, progress, checkpoints)
        continue_solver(t, tf, timer) || break

        evolve_one_time_step!(method, adaptive, controller,
                              FE, y, t, dt, ode_wrap!, update_cache,
                              linear_cache, stage_finder)
        t[1] += dt[1]
        timer.total_steps[1] += 1
    end
    compute_stats!(sol, save_solution, adaptive, timer,
                   stage_finder, loop_stats, config_bytes)
    return nothing
end

"""
    evolve_ode(y0::Union{T, Vector{T}}, t0::T1, tf::Float64, dt0::Float64,
               dy_dt!::Function, options::SolverOptions;
               model_parameters = nothing,
               precision::Type{T2} = Float64) where {T <: AbstractFloat,
                                                     T1 <: AbstractFloat,
                                                     T2 <: AbstractFloat}

Required parameters: `y0`, `dt0`, `dy_dt!`, `options`
"""
function evolve_ode(y0::Union{T, Vector{T}}, t0::T1, tf::Float64, dt0::Float64,
                    dy_dt!::Function, options::SolverOptions;
                    model_parameters = nothing,
                    precision::Type{T2} = Float64) where {T <: AbstractFloat,
                                                          T1 <: AbstractFloat,
                                                          T2 <: AbstractFloat}
    sol = Solution(; precision)
    evolve_ode!(sol, y0, t0, tf, dt0, dy_dt!, options; model_parameters)
    return sol
end
