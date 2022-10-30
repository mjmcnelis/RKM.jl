
# TODO: move to file
function evolve_one_time_step!(method::RungeKutta, adaptive::Fixed,
                               y::Vector{<:AbstractFloat}, t::AbstractFloat, 
                               dt::AbstractFloat, dy_dt!::Function, 
                               dy, y_tmp, f_tmp,
                               )
    butcher = method.butcher
    nrow, ncol = size(butcher) 
    stages = nrow - 1

    # dy = zeros(stages, dimensions)        # TODO: preallocate
    # y_tmp = zeros(dimensions)        
    # f_tmp = zeros(dimensions) 

    c = view(butcher, 1:nrow-1, 1)          # c_i
    A = view(butcher, 1:nrow-1, 2:ncol)     # A_ij
    b = view(butcher, nrow,     2:ncol)     # b_i
    
    # able to kill allocations by making butcher a concrete typed matrix
    for i = 1:stages                        # note: assumes explicit
        t_tmp = t + c[i]*dt
        y_tmp .= y
        for j = 1:i-1
            y_tmp .+= A[i,j] * dy[j]
        end
        dy_dt!(f_tmp, t_tmp, y_tmp)

        dy[i, :] .= dt .* f_tmp
    end    

    for j = 1:stages 
        y .+= b[j] * dy[j]
    end
    nothing
end

function evolve_ode(y0::Union{AbstractFloat, Vector{<:AbstractFloat}}, 
                    t_span::TimeSpan, dy_dt!::Function; 
                    parameters::Parameters, wtime_min = 60)

    @unpack t0, tf, dt0 = t_span
    @unpack adaptive, method = parameters
    
    precision = get_precision(method)

    # initial conditions
    y  = precision[copy(y0)...]
    t  = t0 |> precision
    dt = dt0 |> precision

    # initialize arrays 
    y_arr  = Vector{Vector{precision}}()
    t_arr  = Vector{precision}()
    dt_arr = Vector{precision}()

    time_limit = Dates.now() + Dates.Minute(round(wtime_min))

    # TEMP 
    dimensions = size(y, 1)
    stages = size(method.butcher, 1) - 1
    dy = zeros(stages, dimensions) 
    y_tmp = zeros(dimensions)
    f_tmp = zeros(dimensions)

    # TODO: wrap inside function 
    while true
        push!(y_arr, copy(y)) 
        push!(t_arr, t)
       
        evolve_one_time_step!(method, adaptive, y, t, dt, dy_dt!, dy, y_tmp, f_tmp)

        time_now = Dates.now()

        # TODO: throw LongSolve exception
        (t < tf && time_now < time_limit) || break

        push!(dt_arr, dt)
        t += dt
    end
    
    Solution(; y = y_arr, t = t_arr, dt = dt_arr)
end
