
# TODO: move to file
function evolve_one_time_step!(method::RungeKutta, adaptive::Fixed,
                               y, t, dt, dy_dt!, dy, y_tmp, f_tmp)            
    butcher = method.butcher
    nrow, ncol = size(butcher) 
    stages = nrow - 1

    c = view(butcher, 1:nrow-1, 1)          # c_i
    A = view(butcher, 1:nrow-1, 2:ncol)     # A_ij
    b = view(butcher, nrow,     2:ncol)     # b_i
    
    for i = 1:stages                        # note: assumes explicit
        t_tmp = t + c[i]*dt
        y_tmp .= y
        for j = 1:i-1
            y_tmp .+= A[i,j] * dy[j]
        end
        dy_dt!(f_tmp, t_tmp, y_tmp)

        dy[i,:] .= dt .* f_tmp
    end    

    for j = 1:stages 
        y .+= b[j] * dy[j]
    end
    nothing
end

function evolve_ode(y0, dy_dt!::Function; parameters::Parameters, wtime_min = 1)

    @unpack adaptive, method, t_span = parameters
    @unpack t0, tf, dt0 = t_span
    
    precision = get_precision(method)

    # initial conditions
    y  = precision[copy(y0)...]
    t  = t0  |> precision
    dt = dt0 |> precision
    tf = tf  |> precision

    time_limit = Dates.now() + Dates.Minute(round(wtime_min))

    dimensions = size(y, 1)
    stages = size(method.butcher, 1) - 1
    
    dy    = zeros(stages, dimensions) 
    y_tmp = zeros(dimensions)
    f_tmp = zeros(dimensions)

    # initalize solution
    sol = Solution(; precision) 

    while true
        push!(sol.y, copy(y))
        push!(sol.t, t)
    
        evolve_one_time_step!(method, adaptive, y, t, dt, dy_dt!, dy, y_tmp, f_tmp)

        # TODO: throw LongSolve exception
        (t < tf && Dates.now() < time_limit) || break

        # why is this allocating?
        t += dt
    end
    # TODO: compute dt_arr from t_arr
    sol
end
