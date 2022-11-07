
function evolve_ode(y0, dy_dt!::Function; parameters::Parameters, wtime_min = 1)

    @unpack adaptive, method, t_span = parameters
    @unpack t0, tf, dt0 = t_span
    
    precision = get_precision(method)

    # initial conditions
    y  = precision[copy(y0)...]
    t  = t0  |> precision
    dt = dt0 |> precision
    tf = tf  |> precision
    
    # TODO: may make t_tmp ~ 2D vector
    dt_next = dt

    time_limit = Dates.now() + Dates.Minute(round(wtime_min))

    dimensions = size(y, 1)
    stages = size(method.butcher, 1) - 1
    
    dy    = zeros(stages, dimensions) 
    y_tmp = zeros(dimensions)
    f_tmp = zeros(dimensions)
    f     = zeros(dimensions)

    # TEMP for step doubling
    y1 = zeros(dimensions)
    y2 = zeros(dimensions)
    error = zeros(dimensions)

    # initalize solution
    sol = Solution(; precision) 

    while true
        # TODO: is there a way I can avoid copy(y) without bugs?
        push!(sol.y, copy(y))
        push!(sol.t, t)

        # TODO: preallocate ~ dt_tmp as 2D vector (holds current, and projected time step)
        dt, dt_next = evolve_one_time_step!(method, adaptive, y, t, dt_next, dy_dt!, 
                                            dy, y_tmp, f_tmp, f, y1, y2, error)

        # TODO: split up into two break lines so can throw LongSolve exception
        (t < tf && Dates.now() < time_limit) || break

        # why is this allocating?..
        t += dt
    end
    sol
end
