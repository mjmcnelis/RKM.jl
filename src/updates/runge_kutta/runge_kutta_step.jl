# benchmark.jl note: don't see as much benefit to @muladd as @..
# benchmark.jl note: @.. doesn't help much when switch to MVector
@muladd function runge_kutta_step!(method::RungeKutta, ::Explicit, y::VectorMVector, t::T,
                     dt::T, ode_wrap!::ODEWrapper, dy::MatrixMMatrix, y_tmp::VectorMVector,
                     f_tmp::VectorMVector, FE::MVector{1,Int64},
                     args...) where T <: AbstractFloat
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
        ode_wrap!(f_tmp, t_tmp, y_tmp)
        FE[1] += 1
        @.. dy[:,i] = dt * f_tmp
    end
    @.. y_tmp = y                                        # evaluate iteration
    for j = 1:stages
        dy_stage = view(dy,:,j)
        @.. y_tmp = y_tmp + b[j]*dy_stage
    end
    return nothing
end

@muladd function runge_kutta_step!(method::RungeKutta, ::DiagonalImplicit,
                     y::VectorMVector, t::T, dt::T, ode_wrap!::ODEWrapper,
                     dy::MatrixMMatrix, y_tmp::VectorMVector, f_tmp::VectorMVector,
                     FE::MVector{1,Int64}, error::VectorMVector,
                     J::MatrixMMatrix, linear_cache,
                     stage_finder::ImplicitStageFinder) where T <: AbstractFloat

    @unpack c, A_T, b, stages, explicit_stage, fsal = method
    @unpack root_method, jacobian_method, epsilon, max_iterations, p_norm = stage_finder

    for i = 1:stages
        # first explicit stage should already be pre-evaluated elsewhere
        if i == 1 && explicit_stage[i]
            continue
        end

        # set intermediate time in wrapper
        t_tmp = t + c[i]*dt
        ode_wrap!.t[1] = t_tmp

        # sum over known stages
        @.. y_tmp = y
        for j = 1:i-1
            dy_stage = view(dy,:,j)
            @.. y_tmp = y_tmp + A_T[j,i]*dy_stage
        end

        # guess stage before iterating
        # TODO: look into predictors
        ode_wrap!(f_tmp, t_tmp, y_tmp)
        FE[1] += 1
        @.. dy[:,i] = dt * f_tmp

        if !explicit_stage[i]
            for n = 1:max_iterations
                # reset y_tmp
                @.. y_tmp = y
                for j = 1:i-1
                    dy_stage = view(dy,:,j)
                    @.. y_tmp = y_tmp + A_T[j,i]*dy_stage
                end

                # evaluate current state and ODE
                dy_stage = view(dy,:,i)
                @.. y_tmp = y_tmp + A_T[i,i]*dy_stage
                ode_wrap!(f_tmp, t_tmp, y_tmp)
                FE[1] += 1

                # compute residual error of root equation
                # dy - dt.f(t_tmp, y_tmp + A.dy) = 0
                @.. error = dy_stage - dt*f_tmp

                # compute norms and tolerance
                e_norm  = norm(error, p_norm)           # compute norms
                dy_norm = norm(view(dy,:,i), p_norm)
                tol = epsilon * dy_norm

                # check for root convergence
                # TODO: test doing at least one iteration on helping stiff problems
                if e_norm < tol
                    break
                elseif n == max_iterations - 1
                    # println("failed to converge")
                    # note: allow f_tmp to store dy_dt!
                    # of last iteration for FSAL methods
                    # TODO: count convergence failures in stats
                    break
                end

                if root_method isa FixedPoint
                    @.. dy[:,i] -= error
                elseif root_method isa Newton
                    # evaluate current Jacobian
                    evaluate_system_jacobian!(jacobian_method, FE, J,
                                              ode_wrap!, y_tmp, f_tmp)
                    J .*= (-A_T[i,i]*dt)                 # J <- I - A.dt.J
                    for i in diagind(J)
                        J[i] += 1.0
                    end

                    # pass Jacobian and residual error to linear cache
                    linear_cache.A = J
                    linear_cache.b = error
                    # note: may not need this if use regular newton method
                    # linear_cache = solve_linear_tmp(linear_cache)
                    # @.. dy[:,i] -= linear_cache.u
                    sol = solve(linear_cache)
                    @.. dy[:,i] -= sol.u
                end
            end
        end
    end
    # evaluate update
    @.. y_tmp = y
    for j = 1:stages
        dy_stage = view(dy,:,j)
        @.. y_tmp = y_tmp + b[j]*dy_stage
    end
    return nothing
end