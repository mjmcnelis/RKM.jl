
# TODO: routine only works for explicit
function fixed_runge_kutta_step!(method::RungeKutta, y, t, dt, dy_dt!, dy, y_tmp, f_tmp)
           
    butcher = method.butcher                            # get butcher tableau
    nrow, ncol = size(butcher) 
    stages = nrow - 1                                   # number of stages 

    c = view(butcher, 1:nrow-1, 1)                      # c_i
    A = view(butcher, 1:nrow-1, 2:ncol)                 # A_ij
    b = view(butcher, nrow, 2:ncol)                     # b_i

    for i = 2:stages                                    # evaluate remaining stages
        t_tmp = t + c[i]*dt                             # assumes first stage pre-evaluated
        y_tmp .= y
        for j = 1:i-1
            y_tmp .+= A[i,j] * dy[j]
        end
        dy_dt!(f_tmp, t_tmp, y_tmp)
        dy[i,:] .= dt .* f_tmp
    end    

    y_tmp .= y                                          # evaluate iteration
    for j = 1:stages 
        for i in eachindex(y_tmp)
            y_tmp[i] += b[j] * dy[j,i]
        end
    end
    nothing
end

function evolve_one_time_step!(method::RungeKutta, adaptive::Fixed,
                               y, t, dt, dy_dt!, dy, y_tmp, f_tmp, 
                               f, args...) 
    t = t[1]
    dt = dt[1]
    # TODO: not sure why putting dy_dt! here this kills allocations
    dy_dt!(f, t, y)                                     # evalute first state at (t,y
    dy[1,:] .= dt .* f

    fixed_runge_kutta_step!(method, y, t, dt, dy_dt!, dy, y_tmp, f_tmp)
    y .= y_tmp                                          # get iteration
    nothing
end

function evolve_one_time_step!(method::RungeKutta, adaptive::Doubling,
                               y, t, dt, dy_dt!, dy, y_tmp, f_tmp, 
                               f, y1, y2, error)   
    
    @unpack epsilon, low, high, safety, p_norm, dt_min, dt_max, max_attempts = adaptive

    order = method.order[1]                             # order of scheme
    high ^= order / (1.0 + order)                       # rescale high based on order

    t = t[1]
    dy_dt!(f, t, y)                                     # evaluate first stage at (t,y)

    dt[1] = dt[2]                                       # initialize time step
    rescale = 1.0 

    a = 1
    while true                                          # start step doubling routine 
        # increase dt for next attempt (a > 1)
        dt[1] = min(dt_max, max(dt_min, dt[1]*rescale))       

        dy[1,:] .= dt[1] .* f                          # iterate full time step 
        fixed_runge_kutta_step!(method, y, t, dt[1], dy_dt!, dy, y_tmp, f_tmp)
        y1 .= y_tmp                                     # y1(t+dt)
        
        dy[1,:] .= (dt[1]/2.0) .* f                     # iterate two half time steps
        fixed_runge_kutta_step!(method, y, t, dt[1]/2.0, dy_dt!, dy, y_tmp, f_tmp)
        y2 .= y_tmp                                     # y2(t+dt/2)
        dy_dt!(f_tmp, t + dt[1]/2.0, y2)
        dy[1,:] .= (dt[1]/2.0) .* f_tmp
        fixed_runge_kutta_step!(method, y2, t + dt[1]/2.0, dt[1]/2.0, 
                                dy_dt!, dy, y_tmp, f_tmp)
        y2 .= y_tmp                                     # y2(t+dt)

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