
function evolve_ode(y0::Union{AbstractFloat, Vector{<:AbstractFloat}}, t_span, y_prime; 
                    parameters::Parameters)

    @unpack adaptive, method = parameters
    @unpack t0, tf, dt0 = t_span

    precision = get_precision(method)

    # initial conditions
    y  = precision[copy(y0)...]
    t  = t0 |> precision
    dt = dt0 |> precision

    # initialize arrays 
    y_array  = Vector{Vector{precision}}()
    t_array  = Vector{precision}()
    dt_array = Vector{precision}()

    time_limit = Dates.now() + Dates.Minute(1)

    # TODO: wrap inside function 
    while true
        push!(y_array, y)
        push!(t_array, t)

        # TODO fill out step 
        # y += 1

        time_now = Dates.now()

        # TODO: throw LongSolve exception
        (t < tf && time_now < time_limit) || break

        push!(dt_array, dt)
        t += dt
    end
    
    Solution(; y = y_array, t = t_array, dt = dt_array)
end
