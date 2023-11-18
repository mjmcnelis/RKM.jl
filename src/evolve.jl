# TODO update docstring
"""
    evolve_ode!(sol::Solution, y0::Union{T, Vector{T}}, t0::T1, tf::Float64,
                dt0::Float64, dy_dt!::Function, parameters::Parameters;
                model_parameters = nothing) where {T <: AbstractFloat,
                                                   T1 <: AbstractFloat}

Required parameters: `sol`, `y0`, `dt0`, `dy_dt!`, `parameters`
"""
function evolve_ode!(sol::Solution, y0::Union{T, Vector{T}}, t0::T1, tf::Float64,
                     dt0::Float64, dy_dt!::Function, parameters::Parameters;
                     model_parameters = nothing) where {T <: AbstractFloat,
                                                        T1 <: AbstractFloat}
    config_bytes = @allocated begin
        clear_solution!(sol)
        @unpack precision, FE#=, JE=# = sol

        @unpack adaptive, controller, method, t_range, timer,
                stage_finder, static_array, show_progress, save_solution = parameters

        @set! t_range.t0 = t0
        @set! t_range.tf = tf
        @unpack t0, tf = t_range

        reset_timer!(timer)

        method = reconstruct_method(method, precision)
        if !(adaptive isa Fixed)
            controller = reconstruct_controller(controller, method, adaptive, precision)
        end

        @unpack stages = method

        dimensions = size(y0, 1)
        sol.dimensions .= dimensions

        # initial conditions
        y  = y0 .|> precision
        y isa Vector ? nothing : y = [y]
        static_array ? y = MVector{dimensions}(y...) : nothing

        t0 = rationalize(t0) |> precision
        tf = rationalize(tf) |> precision
        dt0 = rationalize(dt0) |> precision

        t  = precision == BigFloat ? [t0] : MVector{1}(t0)
        dt = precision == BigFloat ? [dt0, dt0] : MVector{2}(dt0, dt0)

        # create ODE wrapper function
        # note: used copy(t) to prevent bug in time accumulation
        ode_wrap! = ODEWrapper(deepcopy(t), model_parameters, dy_dt!)

        # note: should not be SA in general but still may want option if size small
        # note: keep in mind of ForwardDiff issues we had with PaT
        dy    = calloc_matrix(static_array, precision, dimensions, stages)
        y_tmp = calloc_vector(static_array, precision, dimensions)
        f_tmp = calloc_vector(static_array, precision, dimensions)
        f     = calloc_vector(static_array, precision, dimensions)
        # TODO: should have option to make it sparse
        J     = calloc_matrix(static_array, precision, dimensions, dimensions)

        # TEMP for step doubling (embedded too probably)
        y1 = calloc_vector(static_array, precision, dimensions)
        y2 = calloc_vector(static_array, precision, dimensions)
        error = calloc_vector(static_array, precision, dimensions)

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

        # @unpack evaluations = stage_finder.jacobian_method
    end

    if save_solution && adaptive isa Fixed
        sizehint_solution!(sol, t_range, dt0, dimensions)
    end

    loop_stats = @timed while true
        if save_solution
            append!(sol.y, y)
            append!(sol.t, t[1])
        end

        show_progress && monitor_progess(t, progress, checkpoints)
        continue_solver(t, tf, timer) || break

        # TODO: see if can pass kwargs
        evolve_one_time_step!(method, adaptive, controller,
                              FE, y, t, dt, ode_wrap!, dy, y_tmp, f_tmp, f, y1, y2,
                              error, J, linear_cache, stage_finder
                            )
        t[1] += dt[1]

        # note: still get excess allocations (could pass it in above)
        # JE[1] = evaluations[1]
        # q()
    end

    # TODO: look into @code_warntype
    # TODO: wrap into compute_stats function
    # compute_step_rejection_rate!(sol, method, adaptive, timer)
    sol.JE .= stage_finder.jacobian_method.evaluations
    sol.runtime .= loop_stats.time

    solution_bytes = sizeof(sol.y) + sizeof(sol.t)
    sol.solution_size .= format_bytes(solution_bytes)

    if (save_solution && adaptive isa Fixed) || loop_stats.bytes == 0
        sol.excess_memory .= format_bytes(loop_stats.bytes)
    else
        sol.excess_memory .= format_bytes(loop_stats.bytes - solution_bytes)
    end

    sol.config_memory .= format_bytes(config_bytes)
    return nothing
end

"""
    evolve_ode(y0::Union{T, Vector{T}}, t0::T1, tf::Float64, dt0::Float64,
               dy_dt!::Function, parameters::Parameters;
               model_parameters = nothing,
               precision::Type{T2} = Float64) where {T <: AbstractFloat,
                                                     T1 <: AbstractFloat,
                                                     T2 <: AbstractFloat}

Required parameters: `y0`, `dt0`, `dy_dt!`, `parameters`
"""
function evolve_ode(y0::Union{T, Vector{T}}, t0::T1, tf::Float64, dt0::Float64,
                    dy_dt!::Function, parameters::Parameters;
                    model_parameters = nothing,
                    precision::Type{T2} = Float64) where {T <: AbstractFloat,
                                                          T1 <: AbstractFloat,
                                                          T2 <: AbstractFloat}
    sol = Solution(; precision)
    evolve_ode!(sol, y0, t0, tf, dt0, dy_dt!, parameters; model_parameters)
    return sol
end
