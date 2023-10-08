# TODO update docstring
"""
    evolve_ode!(sol::Solution, y0::Union{T, Vector{T}}, dy_dt!::Function;
                model_parameters = nothing,
                static_array::Bool = false, show_progress::Bool = true,
                parameters::Parameters) where {T <: AbstractFloat}

Required parameters: `sol`, `y0`, `dy_dt!`, `parameters`
"""
function evolve_ode!(sol::Solution, y0::Union{T, Vector{T}}, dy_dt!::Function,
                     parameters::Parameters; model_parameters = nothing,
                     static_array::Bool = false,
                     show_progress::Bool = true) where {T <: AbstractFloat}

    clear_solution!(sol)
    @unpack precision, FE#=, JE=#, save_solution = sol

    @unpack adaptive, controller, method, t_range, timer, stage_finder = parameters
    @unpack t0, tf, dt0 = t_range

    # note: if code errors out from bug and doesn't not update
    #       solver_finished variable, then can print time's up warning
    reset_timer!(timer)
    start_timer!(timer)

    method = reconstruct_method(method, precision)
    if !(adaptive isa Fixed)
        controller = reconstruct_controller(controller, method, adaptive, precision)
    end

    @unpack stages, iteration = method

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
    ode_wrap! = ODEWrapper(copy(t), model_parameters, dy_dt!)

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
    linear_cache = set_cacheval(linear_cache, fact)
#     ┌ Warning: set_cacheval is deprecated for mutation on the cache. Use `cache.cacheval = cacheval; cache.isfresh = false
# └ @ LinearSolve ~/.julia/packages/LinearSolve/Rq2MY/src/deprecated.jl:72

    stage_finder = set_jacobian_cache(stage_finder, ode_wrap!, f_tmp, y)

    adaptive isa Fixed ? sizehint_solution!(sol, t_range, dimensions) : nothing

    # for progress meter
    checkpoints = collect(LinRange(t0, tf, 101))[2:end]
    progress = Progress(100)

    # @unpack evaluations = stage_finder.jacobian_method

    stats = @timed allocs = @allocated while true
        if save_solution
            append!(sol.y, y)
            append!(sol.t, t[1])
        end

        show_progress && monitor_progess(t, progress, checkpoints)
        continue_solver(t, tf, timer) || break

        # TODO: see if can pass kwargs
        evolve_one_time_step!(method, iteration, adaptive, controller,
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
    compute_step_rejection_rate!(sol, method, adaptive, timer)
    sol.JE .= stage_finder.jacobian_method.evaluations
    sol.excess_memory .= format_bytes(stats.bytes)
    sol.excess_allocations .= allocs
    sol.runtime .= stats.time
    return nothing
end

"""
    evolve_ode(y0::Union{T, Vector{T}}, dy_dt!::Function;
               static_array::Bool = false, show_progress::Bool = true,
               precision::Type{T2} = Float64, parameters::Parameters,
               jacobian! = jacobian_error) where {T <: AbstractFloat,
                                                  T2 <: AbstractFloat}

Required parameters: `y0`, `dy_dt!`, `parameters`
"""
function evolve_ode(y0::Union{T, Vector{T}}, dy_dt!::Function, parameters::Parameters;
                    model_parameters = nothing, save_solution::Bool = true,
                    static_array::Bool = false, show_progress::Bool = true,
                    precision::Type{T2} = Float64) where {T <: AbstractFloat,
                                                          T2 <: AbstractFloat}

    sol = Solution(; precision, save_solution)
    evolve_ode!(sol, y0, dy_dt!, parameters; model_parameters, static_array,
                                             show_progress)
    return sol
end