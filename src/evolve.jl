
function evolve_ode(y0, dy_dt!::Function; parameters::Parameters, wtime_min::Int64 = 1)

    @unpack adaptive, method, t_span = parameters
    @unpack t0, tf, dt0 = t_span
    
    @unpack precision, iteration = method 

    # initial conditions
    y  = precision[copy(y0)...]
    t  = [t0]
    dt = [dt0, dt0]

    time_limit = Dates.now() + Dates.Minute(round(wtime_min))

    dimensions = size(y, 1)
    stages = size(method.butcher, 2) - 1
    
    dy    = zeros(precision, stages, dimensions) 
    y_tmp = zeros(precision, dimensions)
    f_tmp = zeros(precision, dimensions)
    f     = zeros(precision, dimensions)
   
    # TEMP for step doubling (embedded too probably)
    y1 = zeros(precision, dimensions)
    y2 = zeros(precision, dimensions)
    error = zeros(precision, dimensions)

    # initalize solution
    sol = Solution(; precision) 

    while true
        push!(sol.y, copy(y))
        append!(sol.t, t)

        evolve_one_time_step!(method, iteration, adaptive, y, t, dt, dy_dt!, 
                              dy, y_tmp, f_tmp, f, y1, y2, error)

        check_time(t, tf, time_limit) || break
        t .+= dt[1]
    end
    sol
end

function check_time(t::Vector{Float64}, tf::Float64, time_limit::Dates.DateTime)
    # TODO: split up into two break lines so can throw LongSolve exception
    t[1] < tf && Dates.now() < time_limit
end