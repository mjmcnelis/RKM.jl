# benchmark.jl note: don't see as much benefit to @muladd as @..
# benchmark.jl note: @.. doesn't help much when switch to MVector
@muladd function runge_kutta_step!(method::RungeKutta, ::Explicit, y::VectorMVector, t::T, 
                     dt::T, dy_dt!::F, dy::MatrixMMatrix, y_tmp::VectorMVector, 
                     f_tmp::VectorMVector, FE::MVector{1,Int64}, 
                     args...) where {T <: AbstractFloat, F <: Function}
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
        FE[1] += 1
    end
    @.. y_tmp = y                                        # evaluate iteration
    for j = 1:stages
        dy_stage = view(dy,:,j)
        @.. y_tmp = y_tmp + b[j]*dy_stage
    end
    return nothing
end

@muladd function runge_kutta_step!(method::RungeKutta, ::DiagonalImplicit,
                     y::VectorMVector, t::T, dt::T, dy_dt!::F, dy::MatrixMMatrix,
                     y_tmp::VectorMVector, f_tmp::VectorMVector, FE::MVector{1,Int64},
                     J::MatrixMMatrix, linear_cache, dy_dt_wrap!::ODEWrapper,
                     stage_finder::ImplicitStageFinder) where {T <: AbstractFloat, 
                                                               F <: Function}
    @unpack c, A_T, b, stages, explicit_stage = method
    @unpack root_method, jacobian_method, epsilon, max_iterations = stage_finder

    for i = 1:stages        
        t_tmp = t + c[i]*dt
        dy_dt_wrap!.t[1] = t_tmp                         # set intermediate time in wrapper

        @.. y_tmp = y                                    # sum over known stages
        for j = 1:i-1
            dy_stage = view(dy,:,j)
            @.. y_tmp = y_tmp + A_T[j,i]*dy_stage
        end
        
        if explicit_stage[i]
            dy_dt!(f_tmp, t_tmp, y_tmp)
            @.. dy[:,i] = dt * f_tmp
            FE[1] += 1
        else
            # TODO: look into predictors
            dy_dt!(f_tmp, t_tmp, y_tmp)                  # guess stage before iterating
            @.. dy[:,i] = dt * f_tmp

            for n = 1:max_iterations
                dy_stage = view(dy,:,i)
                @.. y_tmp = y_tmp + A_T[i,i]*dy_stage

                if root_method isa Newton                # evaluate current Jacobian
                    evaluate_system_jacobian!(jacobian_method, J, dy_dt_wrap!, 
                                              y_tmp, f_tmp)
                    J .*= (-A_T[i,i]*dt)                 # J <- I - A.dt.J
                    for i in diagind(J)
                        J[i] += 1.0
                    end
                    linear_cache = set_A(linear_cache, J)# pass Jacobian to linear cache
                end

                dy_dt!(f_tmp, t_tmp, y_tmp)              # evaluate current slope
                FE[1] += 1
                
                @.. y_tmp = y_tmp - A_T[i,i]*dy_stage    # undo addition to y_tmp

                if root_method isa FixedPoint
                    @.. dy[:,i] = dt * f_tmp
                elseif root_method isa Newton
                    # TODO: initialize predictor for dy other than zero,
                    #       and sort out @.. equivalent
                    for k in eachindex(f_tmp)
                        f_tmp[k] = dy[k,i] - dt*f_tmp[k]
                    end

                    # compute residual error of root equation:
                    # dy - dt.f(t_tmp, y_tmp + A.dy) = 0
                    res = norm(f_tmp)

                    dy_norm = norm(view(dy,:,i))         # compute error tolerance 
                    tol = epsilon * dy_norm
               
                    if n > 1 && res < tol                # check if Newton method covnerges
                        # print(n)
                        break
                    end
                    # if n == max_iterations
                    #     @warn "exceeded max Newton iterations at t = $t"
                    #     @show res tol
                    #     println("")
                    # end

                    linear_cache = set_b(linear_cache, f_tmp)
                    # note: may not need this if use regular newton method 
                    linear_cache = solve_linear_tmp(linear_cache)
                    @.. dy[:,i] -= linear_cache.u
                    # sol = solve(linear_cache)
                    # @.. dy[:,i] -= sol.u
                end
            end
        end
    end
    @.. y_tmp = y                                        # evaluate iteration
    for j = 1:stages
        dy_stage = view(dy,:,j)
        @.. y_tmp = y_tmp + b[j]*dy_stage
    end
    nothing
end