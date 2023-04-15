# TODO update docstring
"""
    evolve_ode(y0::Union{T, Vector{T}}, dy_dt!::Function; 
               static_array::Bool = false, show_progress::Bool = true, 
               precision::Type{T2} = Float64, parameters::Parameters, 
               jacobian! = jacobian_error) where {T <: AbstractFloat,
                                                  T2 <: AbstractFloat}

The ODE solver loop. 

Required parameters: `y0`, `dy_dt!`, `parameters`

Note: `jacobian!` argument is temporary
"""
function evolve_ode(y0::Union{T, Vector{T}}, dy_dt!::Function; 
                    static_array::Bool = false, show_progress::Bool = true, 
                    precision::Type{T2} = Float64, parameters::Parameters, 
                    jacobian! = jacobian_error) where {T <: AbstractFloat,
                                                       T2 <: AbstractFloat}

    @unpack adaptive, controller, method, t_range, timer = parameters
    @unpack t0, tf, dt0 = t_range

    # note: if code errors out from bug and doesn't not update
    #       solver_finished variable, then can print time's up warning
    reset_timer!(timer)
    start_timer!(timer)

    method = reconstruct_method(method, precision)
    controller = reconstruct_controller(controller, precision)

    @unpack stages, iteration = method

    dimensions = size(y0, 1)

    # initial conditions
    y  = y0 .|> precision
    y isa Vector ? nothing : y = [y]
    static_array ? y = MVector{dimensions}(y...) : nothing 

    t0 = rationalize(t0) |> precision
    tf = rationalize(tf) |> precision
    dt0 = rationalize(dt0) |> precision

    t  = precision == BigFloat ? [t0] : MVector{1}(t0)
    dt = precision == BigFloat ? [dt0, dt0] : MVector{2}(dt0, dt0)

    # note: should not be SA in general but still may want option if size small
    # note: keep in mind of ForwardDiff issues we had with PaT
    dy    = calloc_matrix(static_array, precision, dimensions, stages)
    y_tmp = calloc_vector(static_array, precision, dimensions)
    f_tmp = calloc_vector(static_array, precision, dimensions)
    f     = calloc_vector(static_array, precision, dimensions)

    # TEMP for step doubling (embedded too probably)
    y1 = calloc_vector(static_array, precision, dimensions)
    y2 = calloc_vector(static_array, precision, dimensions)
    error = calloc_vector(static_array, precision, dimensions)

    # initalize solution
    sol = Solution(; precision, dimensions)
    @unpack FE = sol
   
    # TODO: is resize and indexing faster than append? 
    adaptive isa Fixed ? sizehint_solution!(sol, t_range, dimensions) : nothing

    # for progress meter
    checkpoints = collect(LinRange(t0, tf, 101))[2:end]
    progress = Progress(100)

    allocs = 0
    while true
        append!(sol.y, y) 
        append!(sol.t, t[1])

        show_progress && monitor_progess(t, progress, checkpoints)
        continue_solver(t, tf, timer) || break

        # TODO: see if can pass kwargs
        allocs += @allocated evolve_one_time_step!(method, iteration, adaptive, controller, FE, y, t, dt, 
                              dy_dt!, dy, y_tmp, f_tmp, f, y1, y2, error, jacobian!)
        t[1] += dt[1]
    end
    @info "Number of memory allocations during solve = $allocs"

    compute_step_rejection_rate!(sol, method, adaptive, timer)
    return sol
end
