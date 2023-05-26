



# TODO: combine explicit and implicit routines
function evolve_one_time_step!(method::RungeKutta, iteration::DiagonalImplicit,
             adaptive::Embedded, controller::Controller, FE::MVector{1,Int64},  
             y::VectorMVector, t::VectorMVector{1,T}, dt::VectorMVector{2,T},
             dy_dt!::F, dy::MatrixMMatrix, y_tmp::VectorMVector, 
             f_tmp::VectorMVector, f::VectorMVector, y1::VectorMVector, 
             y2::VectorMVector, error::VectorMVector, 
             J::MatrixMMatrix, linear_cache, dy_dt_wrap!, 
             stage_finder::ImplicitStageFinder) where {T <: AbstractFloat, F}

    @unpack epsilon, low, high, safety, p_norm, dt_min, dt_max, max_attempts = adaptive

    y_norm = norm(y, p_norm)                            # compute norm of current state

    order_max = maximum(method.order)                   # max/min orders in embedded scheme
    order_min = minimum(method.order)

    # note: comment if benchmark against OrdinaryDiffEq (add boolean?)
    high    ^= (order_min / order_max)                  # rescale high, epsilon parameters
    epsilon ^= (order_min / order_max)

    dt[1] = dt[2]                                       # initialize time step
    rescale = 1.0

    # # TODO: dispatch fsal
    # if !(method.fsal isa FSAL && FE[1] > 0)
    #     dy_dt!(f, t[1], y)                              # evaluate first stage at (t,y)
    # endz

    attempts = 1
    while true                                          # start step doubling routine
        dt[1] = min(dt_max, max(dt_min, dt[1]*rescale)) # increase dt for next attempt
                       
        runge_kutta_step!(method, iteration, y, t[1], dt[1], dy_dt!, dy, y_tmp, 
                          f_tmp, J, linear_cache, dy_dt_wrap!, stage_finder)
        @.. y1 = y_tmp                                  # primary iteration

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

        if FE[1] == 0                                   # initialize controller
            initialize_controller!(controller, e_norm, dt[1])
        end

        if e_norm == 0.0                                # compute scaling factor for dt
            rescale = high
        else
            rescale = rescale_time_step(controller, tol, e_norm, order_min)
            rescale = min(high, max(low, safety*rescale))
        end

        dt[2] = min(dt_max, max(dt_min, dt[1]*rescale)) # projected dt for next iteration

        if e_norm <= tol                                # compare error to tolerance
            set_previous_control_vars!(controller, e_norm, dt[1])
            break 
        end
        attempts <= max_attempts || (@warn "embedded exceeded $max_attempts attempts"; break)
        attempts += 1
    end
    # # TODO: dispatch fsal
    # if method.fsal isa FSAL
    #     @.. f = f_tmp 
    # end
    
    @.. y = y1                                          # get iteration
    # add_function_evaluations!(FE, iteration, adaptive, method, attempts)
    FE[1] += 10
    return nothing
end
