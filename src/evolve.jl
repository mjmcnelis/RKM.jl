
function evolve_ode(y0::Union{AbstractFloat, Vector{<:AbstractFloat}}, t_span, y_prime; 
                    parameters::Parameters)

    @unpack adaptive, method = parameters
    @unpack t0, tf, dt0 = t_span

    precision = get_precision(method)

    # initial conditions
    # y  = copy(y0) |> precision 
    y  = precision[copy(y0)...]
    t  = t0 |> precision
    dt = dt0 |> precision

    # dimension = size(y, 1)

    # initialize arrays 
    # y_array  = Matrix{precision}(undef, 1, dimension)
    # y_array[1, :] .= y
    y_array  = Vector{Vector{precision}}()
    t_array  = Vector{precision}()
    dt_array = Vector{precision}()

    for i = 1:10
    # while true
        # y_array = hcat(y_array, y)  # terrible allocation
        # push!(y_array, [y,])
        # y_array = [y_array; reshape([y,], dimension, 1)]
        push!(y_array, y)             # can do better if use vector x vector?
        push!(t_array, t)

        # TODO fill out step 
        # y += 1

        # y_array = [y_array; y]      # N x dimension array seems to allocate less
        
        push!(dt_array, dt)

        # TODO: make timelimit
        # t <= tf || break

        # t += dt
    end
    
    Solution(; y = y_array, t = t_array, dt = dt_array)
end
