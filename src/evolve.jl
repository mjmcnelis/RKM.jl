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
    y = MVector{dimensions, precision}(copy(y0)...)

    # y  = y0 .|> precision
    # y isa Vector ? nothing : y = [y]
 
    t  = MVector{1,precision}(t0)
    dt = MVector{2,precision}(dt0, dt0)
  
    # note: should not be MVector in general but still may want option if size small
    # note: keep in mind of ForwardDiff issues we had with PaT
    dy = @MMatrix zeros(precision, dimensions, stages) 
    y_tmp = @MVector zeros(precision, dimensions)
    f_tmp = @MVector zeros(precision, dimensions)
    f = @MVector zeros(precision, dimensions)  # want to test it out though
    
    # TEMP for step doubling (embedded too probably)
    y1 = zeros(precision, dimensions)
    y2 = zeros(precision, dimensions)
    error = zeros(precision, dimensions)

    # initalize solution
    sol = Solution(; precision, dimensions)
    @unpack FE = sol

    # TODO: is resize and fill faster than append? 
    adaptive isa Fixed ? sizehint_solution!(sol, t_range, dimensions) : nothing

    while true
        # update_solution!(sol, y, t)
        append!(sol.y, y)
        append!(sol.t, t[1])

        # TODO: see if can pass kwargs
        evolve_one_time_step!(method, iteration, adaptive, FE, y, t, dt, dy_dt!,
                              dy, y_tmp, f_tmp, f, y1, y2, error, jacobian!)

        continue_solver(t, tf, timer) || break
        @.. t += dt[1]

        # propagate_time(t, dt)
    end
    sol
end

function propagate_time(t::MVector{1,T}, dt::MVector{2,T}) where T <: AbstractFloat
    @.. t += dt[1]
    nothing 
end
