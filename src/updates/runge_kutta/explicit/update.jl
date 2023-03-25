# benchmark.jl note: don't see as much benefit to @muladd as @..
# benchmark.jl note: @.. doesn't help much when switch to MVector
@muladd function fixed_runge_kutta_step!(method::RungeKutta, ::Explicit, 
                     y::VectorMVector, t::T, dt::T, dy_dt!::F, dy::MatrixMMatrix, 
                     y_tmp::VectorMVector, f_tmp::VectorMVector) where {T <: AbstractFloat,
                                                                        F <: Function}
    @unpack c, A_T, b, stages = method

    for i = 2:stages                                    # evaluate remaining stages
        t_tmp = t + c[i]*dt                             # assumes first stage pre-evaluated
        @.. y_tmp = y
        # TODO: need a better dy cache for performance
        for j = 1:i-1
            # TODO: continue if A,b = 0?
            dy_stage = view(dy,:,j)
            @.. y_tmp = y_tmp + A_T[j,i]*dy_stage
        end
        # TODO: skip if intermediate update not needed in next row(s)? 
        dy_dt!(f_tmp, t_tmp, y_tmp)
        @.. dy[:,i] = dt * f_tmp
    end

    @.. y_tmp = y                                        # evaluate iteration
    for j = 1:stages
        dy_stage = view(dy,:,j)
        @.. y_tmp = y_tmp + b[j]*dy_stage
    end
    return nothing
end

function evolve_one_time_step!(method::RungeKutta, iteration::Explicit,
             adaptive::Fixed, FE::MVector{1,Int64}, y::VectorMVector, 
             t::VectorMVector{1,T}, dt::VectorMVector{2,T}, dy_dt!::F,
             dy::MatrixMMatrix, y_tmp::VectorMVector, 
             f_tmp::VectorMVector, args...) where {T <: AbstractFloat, F}
             
    dy_dt!(f_tmp, t[1], y)                              # evaluate first stage at (t,y)
    @.. dy[:,1] = dt[1] * f_tmp

    fixed_runge_kutta_step!(method, iteration, y, t[1], dt[1], dy_dt!, dy, y_tmp, f_tmp)
    @.. y = y_tmp                                       # get iteration

    add_function_evaluations!(FE, iteration, adaptive, method)
    return nothing
end

function evolve_one_time_step!(method::RungeKutta, iteration::Explicit,
             adaptive::FiniteDiff, FE::MVector{1,Int64}, y::Vector{T}, 
             t::Union{Vector{T}, MVector{1,T}}, dt::Union{Vector{T}, MVector{2,T}}, 
             dy_dt!::F, dy::Matrix{T}, y_tmp::Vector{T}, f_tmp::Vector{T}, 
             # TODO: may want to do kwargs for different caches used in adaptive methods
             # note: next argument was f but renamed it to y_prev here
             y_prev::Vector{T}, args...) where {T <: AbstractFloat, F}

    dy_dt!(f_tmp, t[1], y)                              # evaluate first stage at (t,y)

    # TEMP: fixed time step for first update until estimate first time step w/ doubling
    if FE[1] > 0
        @unpack epsilon, low, high, p_norm, dt_min, dt_max = adaptive
    
        order = method.order[1]                         # order of scheme
    
        # for high did 1/order, but epsilon 2/x?
        # don't remember my reasoning for that
        high    ^= 1.0 / order                          # rescale high based on order
        epsilon ^= 2.0 / (1.0 + order)

        # TODO: check if this is allocating (don't set new time step)
        # note: seems like it's allocating even for fixed time step 

        # @.. y_tmp = y + dt[1]*f_tmp                     # compute y_star (stored in y_tmp)
    
        # approximate C w/ central differences (stored in y_tmp)
        # @.. y_tmp = y_tmp - 2.0*y + y_prev
        @.. y_tmp = dt[1]*f_tmp - y + y_prev
        @.. y_tmp *= 2.0/dt[1]^2

        # dy_dt!(y_tmp, t[1], y_prev)
        # @.. y_tmp = f_tmp - y_tmp
        # @.. y_tmp *= 2.0/dt[1]
        # C = 2/dt^2 * (dt*f - y + y_prev)
        # C = 2/dt * (f - (y - y_prev)/dt ) ~ 2/dt * (f_n - f_n-1)
      
        C_norm = norm(y_tmp, p_norm)
        y_norm = norm(y, p_norm)
        f_norm = norm(f_tmp, p_norm)

        # TODO: solve more complicated algebraic equation from VAH paper
        if C_norm == 0.0                                # compute new time step in dt[2]
            dt[2] = high*dt[1]
        else
            if C_norm*y_norm > 2.0*epsilon*f_norm^2
                dt[2] = sqrt(2.0*epsilon*y_norm/C_norm)
            else
                dt[2] = 2.0*epsilon*f_norm/C_norm
            end
        end
        dt[2] = min(high*dt[1], max(low*dt[1], dt[2]))  # control growth rate
        dt[2] = min(dt_max, max(dt_min, dt[2]))         # impose min/max bounds
    end
    # evaluate first stage iteration w/ new time step (i.e. dt[2])
    @.. dy[:,1] = dt[2] * f_tmp    
    fixed_runge_kutta_step!(method, iteration, y, t[1], dt[2], dy_dt!, dy, y_tmp, f_tmp)

    dt[1] = dt[2]                                       # store current time step
    @.. y_prev = y                                      # store current solution
    @.. y = y_tmp                                       # get iteration

    add_function_evaluations!(FE, iteration, adaptive, method)
    return nothing 
end

function evolve_one_time_step!(method::RungeKutta, iteration::Explicit,
             adaptive::Doubling, controller::Controller, FE::MVector{1,Int64},  
             y::Vector{T}, t::Union{Vector{T}, MVector{1,T}}, 
             dt::Union{Vector{T}, MVector{2,T}}, dy_dt!::F, dy::Matrix{T}, 
             y_tmp::Vector{T}, f_tmp::Vector{T}, f::Vector{T}, y1::Vector{T}, 
             y2::Vector{T}, error::Vector{T}, args...) where {T <: AbstractFloat, F}

    @unpack epsilon, low, high, safety, p_norm, dt_min, dt_max, max_attempts = adaptive

    order = method.order[1]                             # order of scheme
    high ^= order / (1.0 + order)                       # rescale high based on order
    epsilon ^= order / (1.0 + order)

    dy_dt!(f, t[1], y)                                  # evaluate first stage at (t,y)

    dt[1] = dt[2]                                       # initialize time step
    rescale = 1.0

    attempts = 1
    while true                                          # start step doubling routine
        dt[1] = min(dt_max, max(dt_min, dt[1]*rescale)) # increase dt for next attempt

        doubling_runge_kutta_step!(method, iteration, y, t[1], dt[1],
                                   dy_dt!, dy, y_tmp, f_tmp, f, y1, y2)

        @.. error = (y2 - y1) / (2.0^order - 1.0)       # estimate local truncation error
        @.. y2 = y2 + error                             # Richardson extrapolation

        e_norm = norm(error, p_norm)                    # compute norms
        y_norm = norm(y2, p_norm)
        Δy = y1
        @.. Δy = y2 - y
        Δy_norm = norm(Δy, p_norm)

        tol = epsilon * max(y_norm, Δy_norm)            # compute tolerance

        if FE[1] == 0                                   # initialize controller
            set_previous_control_error!(controller, e_norm)
        end

        if e_norm == 0.0                                # compute scaling factor for dt
            rescale = 1.0
        else
            rescale = rescale_time_step(controller, tol, e_norm, order)
            rescale = min(high, max(low, safety*rescale))
        end

        dt[2] = min(dt_max, max(dt_min, dt[1]*rescale)) # projected dt for next iteration

        if e_norm <= tol                                # compare error to tolerance
            set_previous_control_error!(controller, e_norm)
            break 
        end
        attempts <= max_attempts || (@warn "step doubling exceeded $max_attempts attempts"; break)
        attempts += 1
    end
    @.. y = y2                                          # get iteration
    add_function_evaluations!(FE, iteration, adaptive, method, attempts)
    return nothing
end

function evolve_one_time_step!(method::RungeKutta, iteration::Explicit,
             adaptive::Embedded, FE::MVector{1,Int64}, y::VectorMVector, 
             t::VectorMVector{1,T}, dt::VectorMVector{2,T}, dy_dt!::F,
             dy::MatrixMMatrix, y_tmp::VectorMVector, f_tmp::VectorMVector, 
             f::VectorMVector, y1::VectorMVector, y2::VectorMVector, error::VectorMVector,
             args...) where {T <: AbstractFloat, F}
           
    @unpack epsilon, low, high, safety, p_norm, dt_min, dt_max, max_attempts = adaptive

    y_norm = norm(y, p_norm)                            # compute norm of current state

    order_max = maximum(method.order)                   # max/min orders in embedded scheme
    order_min = minimum(method.order)

    # note: comment if benchmark against OrdinaryDiffEq (add boolean?)
    high    ^= (order_min / order_max)                  # rescale high, epsilon parameters
    epsilon ^= (order_min / order_max)

    # TODO: dispatch fsal
    if !(method.fsal isa FSAL && FE[1] > 0)
        dy_dt!(f, t[1], y)                                  # evaluate first stage at (t,y)
    end

    dt[1] = dt[2]                                       # initialize time step

    # TODO: float type can change from Float64 to Double64
    # note: is that why benchmark vs OrdinaryDiffEq bad for Double64?
    rescale = 1.0

    attempts = 1
    while true                                          # start embedded routine
        dt[1] = min(dt_max, max(dt_min, dt[1]*rescale)) # increase dt for next attempt

        @.. dy[:,1] = dt[1] * f                         # primary iteration
        fixed_runge_kutta_step!(method, iteration, y, t[1], dt[1],
                                dy_dt!, dy, y_tmp, f_tmp)
        @.. y1 = y_tmp

        embedded_runge_kutta_step!(method, y, dy, y_tmp)
        @.. y2 = y_tmp                                  # embedded iteration

        @.. error = y2 - y1                             # local error of embedded pair

        e_norm = norm(error, p_norm)                    # compute norms
        y1_norm = norm(y1, p_norm)
        # TODO: need to use Δy = y2 since it's secondary method
        #       but should make labeling consistent w/ doubling
        Δy = y2
        @.. Δy = y1 - y
        Δy_norm = norm(Δy, p_norm)

        # compute tolerance
        tol = epsilon * max(y_norm, y1_norm, Δy_norm)

        if e_norm == 0.0                                # compute scaling factor for dt
            rescale = high          
        else
            rescale = (tol / e_norm)^(1.0/(1.0 + order_min))
            rescale = min(high, max(low, safety*rescale))
        end

        dt[2] = min(dt_max, max(dt_min, dt[1]*rescale)) # projected dt for next iteration

        e_norm > tol || break                           # compare error to tolerance
        attempts <= max_attempts || (@warn "embedded exceeded $max_attempts attempts"; break)
        attempts += 1
    end
    
    # TODO: dispatch fsal
    if method.fsal isa FSAL
        @.. f = f_tmp 
    end

    @.. y = y1                                          # get iteration
    add_function_evaluations!(FE, iteration, adaptive, method, attempts)
    return nothing
end

function doubling_runge_kutta_step!(method, iteration::Explicit, y, t, dt,
                                    dy_dt!, dy, y_tmp, f_tmp, f, y1, y2)
    @.. dy[:,1] = dt * f                                # iterate full time step
    fixed_runge_kutta_step!(method, iteration, y, t, dt, dy_dt!, dy, y_tmp, f_tmp)
    @.. y1 = y_tmp                                      # y1(t+dt)

    @.. dy[:,1] = (dt/2.0) * f                          # iterate two half time steps
    fixed_runge_kutta_step!(method, iteration, y, t, dt/2.0, dy_dt!, dy, y_tmp, f_tmp)
    @.. y2 = y_tmp                                      # y2(t+dt/2)
    dy_dt!(f_tmp, t + dt/2.0, y2)
    @.. dy[:,1] = (dt/2.0) * f_tmp
    fixed_runge_kutta_step!(method, iteration, y2, t + dt/2.0, dt/2.0,
    dy_dt!, dy, y_tmp, f_tmp)
    @.. y2 = y_tmp                                      # y2(t+dt)
    return nothing
end

@muladd function embedded_runge_kutta_step!(method, y, dy, y_tmp)
    @unpack stages, b_hat = method
    @.. y_tmp = y                                       # evaluate iteration
    for j = 1:stages
        dy_stage = view(dy,:,j)
        @.. y_tmp = y_tmp + b_hat[j]*dy_stage
    end
    return nothing
end
