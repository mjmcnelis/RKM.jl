# TODO update docstring
"""
    evolve_ode(y0::Union{T, Vector{T}}, dy_dt!::Function; jacobian! = jacobian_error,
               parameters::Parameters
               precision::Type{T2} = Float64) where {T <: AbstractFloat,
                                                     T2 <: AbstractFloat}

The ODE solver loop. 

Required parameters: `y0`, `dy_dt!`, `parameters`

Note: `jacobian!` argument is temporary
"""
function evolve_ode(y0::Union{T, Vector{T}}, dy_dt!::Function; jacobian! = jacobian_error,
                    parameters::Parameters,  
                    precision::Type{T2} = Float64) where {T <: AbstractFloat,
                                                          T2 <: AbstractFloat}

    @unpack adaptive, method, t_range, timer = parameters
    @unpack t0, tf, dt0 = t_range

    method = reconstruct_method(method, precision)
    @unpack stages, iteration = method

    dimensions = size(y0, 1)

    # initial conditions
    y  = y0 .|> precision
    y isa Vector ? nothing : y = [y]

    # note: testing 
    t0 = rationalize(t0)
    tf = rationalize(tf) |> precision
    dt0 = rationalize(dt0)

    t  = precision == BigFloat ? precision[t0] : MVector{1,precision}(t0)
    dt = precision == BigFloat ? precision[dt0, dt0] : MVector{2,precision}(dt0, dt0)

    # note: should not be SA in general but still may want option if size small
    # note: keep in mind of ForwardDiff issues we had with PaT
    dy    = zeros(precision, dimensions, stages)
    y_tmp = zeros(precision, dimensions)
    f_tmp = zeros(precision, dimensions)
    f     = zeros(precision, dimensions)

    # TEMP for step doubling (embedded too probably)
    y1 = zeros(precision, dimensions)
    y2 = zeros(precision, dimensions)
    error = zeros(precision, dimensions)

    # initalize solution
    sol = Solution(; precision, dimensions)
    @unpack FE = sol

    # TODO: is resize and indexing faster than append? 
    adaptive isa Fixed ? sizehint_solution!(sol, t_range, dimensions) : nothing

    while true
        append_solution!(sol, y, t)
        
        continue_solver(t, tf, timer) || break

        # TODO: see if can pass kwargs
        evolve_one_time_step!(method, iteration, adaptive, FE, y, t, dt, dy_dt!,
                              dy, y_tmp, f_tmp, f, y1, y2, error, jacobian!)

        @.. t += dt[1]
        # TODO: make a function monitoring evolution progress
    end
    compute_step_rejection_rate!(sol, method, adaptive, timer)
    return sol
end

# TODO: move to utils?
function rationalize(x; sigdigits = 16)
    Int(round(x*10^(sigdigits-1),digits=0))//10^(sigdigits-1)
end
