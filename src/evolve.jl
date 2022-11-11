
function evolve_ode(y0, dy_dt!::Function; parameters::Parameters, wtime_min = 1)

    @unpack adaptive, method, t_span = parameters
    @unpack t0, tf, dt0 = t_span
    
    precision = method.precision 

    # initial conditions
    y  = precision[copy(y0)...]
    t  = t0  |> precision
    dt = dt0 |> precision
    tf = tf  |> precision

    time_limit = Dates.now() + Dates.Minute(round(wtime_min))

    dimensions = size(y, 1)
    stages = size(method.butcher, 1) - 1
    
    dy     = zeros(stages, dimensions) 
    y_tmp  = zeros(dimensions)
    f_tmp  = zeros(dimensions)
    f      = zeros(dimensions)
    t_vec  = zeros(1)
    t_vec .= t
    dt_vec = zeros(2)                       # holds [dt_current, dt_next]
    dt_vec .= dt
    
    # TEMP for step doubling
    y1 = zeros(dimensions)
    y2 = zeros(dimensions)
    error = zeros(dimensions)

    # initalize solution
    sol = Solution(; precision) 

    while true
        push!(sol.y, copy(y))
        append!(sol.t, t_vec)

        evolve_one_time_step!(method, adaptive, y, t_vec, dt_vec, dy_dt!, 
                              dy, y_tmp, f_tmp, f, y1, y2, error)

        # TODO: split up into two break lines so can throw LongSolve exception
        (t_vec[1] < tf && Dates.now() < time_limit) || break

        t_vec .+= dt_vec[1]
    end
    sol
end
