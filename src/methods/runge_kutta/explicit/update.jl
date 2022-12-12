
# TODO: so far, routine only works for an explicit, primary method
function fixed_runge_kutta_step!(method::RungeKutta, ::Explicit, 
             y::Vector{T}, t::Float64, dt::Float64, dy_dt!::Function, 
             dy::Matrix{T}, y_tmp::Vector{T}, f_tmp::Vector{T}) where {T <: AbstractFloat}

    @unpack c, A, b, stages = method
    
    for i = 2:stages                                    # evaluate remaining stages
        t_tmp = t + c[i]*dt                             # assumes first stage pre-evaluated
        y_tmp .= y
        for j = 1:i-1
            y_tmp .+= A[i,j] .* view(dy, j, :)
        end
        dy_dt!(f_tmp, t_tmp, y_tmp)
        dy[i,:] .= dt .* f_tmp 
    end 

    y_tmp .= y                                          # evaluate iteration
    for j = 1:stages 
        y_tmp .+= b[j] .* view(dy, j, :)
    end
    nothing
end

function embedded_runge_kutta_step!(method, y, dy, y_tmp)
    @unpack stages, b_hat = method
    y_tmp .= y                                          # evaluate iteration
    for j = 1:stages 
        y_tmp .+= b_hat[j] .* view(dy, j, :)
    end
    nothing
end

function doubling_runge_kutta_step!(method, iteration::Explicit, y, t, dt,
                                    dy_dt!, dy, y_tmp, f_tmp, f, y1, y2)
    dy[1,:] .= dt .* f                                  # iterate full time step 
    fixed_runge_kutta_step!(method, iteration, y, t, dt, dy_dt!, dy, y_tmp, f_tmp)
    y1 .= y_tmp                                         # y1(t+dt)
    
    dy[1,:] .= (dt/2.0) .* f                            # iterate two half time steps
    fixed_runge_kutta_step!(method, iteration, y, t, dt/2.0, dy_dt!, dy, y_tmp, f_tmp)
    y2 .= y_tmp                                         # y2(t+dt/2)
    dy_dt!(f_tmp, t + dt/2.0, y2)
    dy[1,:] .= (dt/2.0) .* f_tmp
    fixed_runge_kutta_step!(method, iteration, y2, t + dt/2.0, dt/2.0, dy_dt!, dy, y_tmp,
                            f_tmp)
    y2 .= y_tmp                                         # y2(t+dt)
    nothing
end

function evolve_one_time_step!(method::RungeKutta, iteration::Explicit, ::Fixed,
             y::Vector{T}, t::MVector{1,Float64}, dt::MVector{2,Float64}, dy_dt!::Function,
             # TODO: try F concrete type instead of Function
             dy::Matrix{T}, y_tmp::Vector{T}, f_tmp::Vector{T}, f::Vector{T}, 
             args...) where {T <: AbstractFloat}
    # TODO: not sure why putting dy_dt! here this kills allocations
    dy_dt!(f, t[1], y)                                  # evalute first state at (t,y)
    dy[1,:] .= dt[1] .* f

    fixed_runge_kutta_step!(method, iteration, y, t[1], dt[1], dy_dt!, dy, y_tmp, f_tmp)

    y .= y_tmp                                          # get iteration
    nothing
end

function evolve_one_time_step!(method::RungeKutta, iteration::Explicit, adaptive::Doubling,
             y::Vector{T}, t::MVector{1,Float64}, dt::MVector{2,Float64}, dy_dt!::Function,
             dy::Matrix{T}, y_tmp::Vector{T}, f_tmp::Vector{T}, f::Vector{T},
             y1::Vector{T}, y2::Vector{T}, error::Vector{T}, 
             args...) where {T <: AbstractFloat}
    
    @unpack epsilon, low, high, safety, p_norm, dt_min, dt_max, max_attempts = adaptive

    order = method.order[1]                             # order of scheme
    high ^= order / (1.0 + order)                       # rescale high based on order

    dy_dt!(f, t[1], y)                                  # evaluate first stage at (t,y)

    dt[1] = dt[2]                                       # initialize time step
    rescale = 1.0 

    a = 1
    while true                                          # start step doubling routine 
        dt[1] = min(dt_max, max(dt_min, dt[1]*rescale)) # increase dt for next attempt      

        doubling_runge_kutta_step!(method, iteration, y, t[1], dt[1], 
                                   dy_dt!, dy, y_tmp, f_tmp, f, y1, y2)

        error .= (y2 .- y1) ./ (2.0^order - 1.0)        # estimate local truncation error
        y2 .+= error                                    # Richardson extrapolation

        e_norm = LinearAlgebra.norm(error, p_norm)      # compute norms
        y_norm = LinearAlgebra.norm(y2, p_norm)
        Δy = y1
        Δy .= y2 .- y
        Δy_norm = LinearAlgebra.norm(Δy, p_norm)

        tol = epsilon * max(y_norm, Δy_norm)            # compute tolerance

        if e_norm == 0.0                                # compute scaling factor for dt 
            rescale = 1.0 
        else
            rescale = (tol / e_norm)^(1.0/(1.0 + order))
            rescale = min(high, max(low, safety*rescale))
        end
        
        dt[2] = min(dt_max, max(dt_min, dt[1]*rescale)) # projected dt for next iteration
        
        e_norm > tol || break                           # compare error to tolerance
        a <= max_attempts || (@warn "step doubling exceeded $max_attempts attempts"; break)
        a += 1
    end
    y .= y2
    nothing
end

function evolve_one_time_step!(method::RungeKutta, iteration::Explicit, adaptive::Embedded,
             y::Vector{T}, t::MVector{1,Float64}, dt::MVector{2,Float64}, dy_dt!::Function,
             dy::Matrix{T}, y_tmp::Vector{T}, f_tmp::Vector{T}, f::Vector{T},
             y1::Vector{T}, y2::Vector{T}, error::Vector{T}, 
             args...) where {T <: AbstractFloat}

    @unpack epsilon, low, high, safety, p_norm, dt_min, dt_max, max_attempts = adaptive

    order_max = maximum(method.order)                   # max/min orders in embedded scheme
    order_min = minimum(method.order)

    high    ^= (order_min / order_max)                  # rescale high, epsilon parameters
    # this caused issues (look into this)
    epsilon ^= (order_min / order_max)

    dy_dt!(f, t[1], y)                                  # evaluate first stage at (t,y)
    
    dt[1] = dt[2]                                       # initialize time step
    rescale = 1.0

    a = 1
    while true                                          # start embedded routine 
        dt[1] = min(dt_max, max(dt_min, dt[1]*rescale)) # increase dt for next attempt   

        dy[1,:] .= dt[1] .* f                           # primary iteration
        fixed_runge_kutta_step!(method, iteration, y, t[1], dt[1], 
                                dy_dt!, dy, y_tmp, f_tmp)
        y1 .= y_tmp

        embedded_runge_kutta_step!(method, y, dy, y_tmp)
        y2 .= y_tmp                                     # embedded iteration

        error .= (y2 .- y1)                             # local error of embedded pair

        e_norm = LinearAlgebra.norm(error, p_norm)      # compute norms
        y_norm = LinearAlgebra.norm(y1, p_norm)
        # TODO: need to use Δy =  y2 since it's secondary method 
        #       but should make labeling consistent w/ doubling
        Δy = y2                         
        Δy .= y1 .- y
        Δy_norm = LinearAlgebra.norm(Δy, p_norm)

        tol = epsilon * max(y_norm, Δy_norm)            # compute tolerance

        if e_norm == 0.0                                # compute scaling factor for dt
            rescale = 1.0 
        else
            rescale = (tol / e_norm)^(1.0/(1.0 + order_min))
            rescale = min(high, max(low, safety*rescale))
        end
        
        dt[2] = min(dt_max, max(dt_min, dt[1]*rescale)) # projected dt for next iteration
        
        e_norm > tol || break                           # compare error to tolerance
        a <= max_attempts || (@warn "step doubling exceeded $max_attempts attempts"; break)
        a += 1
    end
    y .= y1
    nothing
end
