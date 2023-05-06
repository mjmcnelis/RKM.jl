
@muladd function fixed_runge_kutta_step!(method::RungeKutta, ::DiagonalImplicit,
                     y::VectorMVector, t::T, dt::T, dy_dt!::F, dy::MatrixMMatrix,
                     y_tmp::VectorMVector, f_tmp::VectorMVector, J::MatrixMMatrix,
                     linear_cache, finitediff_cache, jacobian_config, 
                     dy_dt_wrap!) where {T <: AbstractFloat, F <: Function}

    @unpack c, A_T, b, stages = method

    # root_solver = "fixed_point"        # will just use fixed point iteration for now
    root_solver = "newton_fast"
    eps_root = 1e-8
    max_iterations = 2                   # 2 for benchmarking ImplicitEuler(), fixed dt

    for i = 1:stages
        # evaluate jacobian
        if root_solver == "newton_fast" 
            # TODO: how to reduce allocations here, take it apart or make a wrapper?
            # so in order for jacobian config to work, dy_dt_wrap! argument
            # has to be the same object stored in jacobian_config
            jacobian!(J, dy_dt_wrap!, f_tmp, y, jacobian_config)
            # @show J
            # finite_difference_jacobian!(J, dy_dt_wrap!, y, finitediff_cache) 
            # @show J; q()

            J .*= (-A_T[i,i]*dt)
            for i in diagind(J)
                J[i] += 1.0
            end
            # pass jacobian to linear cache
            linear_cache = set_A(linear_cache, J)
        end
        t_tmp = t + c[i]*dt
        # sum over known stages
        @.. y_tmp = y 
        for j = 1:i-1
            dy_stage = view(dy,:,j)
            @.. y_tmp = y_tmp + A_T[j,i]*dy_stage
        end
      
        # TEMP iterate w/o any breaks for now
        for n = 1:max_iterations
            dy_stage = view(dy,:,i)
            @.. y_tmp = y_tmp + A_T[i,i]*dy_stage
    
            # evaluate ODE at previous iteration (could maybe just iterate y_tmp itself)
            dy_dt!(f_tmp, t_tmp, y_tmp)

            # TEMP undo addition (for minimizing allocations)
            @.. y_tmp = y_tmp - A_T[i,i]*dy_stage

            if root_solver == "fixed_point"
                @.. dy[:,i] = dt * f_tmp
            elseif root_solver == "newton_fast"
                # TODO: try to solve for dy directly instead of d(dy)
                for k in eachindex(f_tmp)
                    f_tmp[k] = dy[k,i] - dt*f_tmp[k]
                end
                # from python
                # g = z - dt*y_prime(t + dt*c[i], y + dy + z*Aii)

                linear_cache = set_b(linear_cache, f_tmp)
                linear_cache = solve_linear_tmp(linear_cache)
                @.. dy[:,i] -= linear_cache.u
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

# TODO: use alternative indents
function evolve_one_time_step!(method::RungeKutta, iteration::DiagonalImplicit,
            adaptive::Fixed, ::Controller, FE::MVector{1,Int64},
            y::VectorMVector, t::VectorMVector{1,T}, dt::VectorMVector{2,T},
            dy_dt!::F, dy::MatrixMMatrix, y_tmp::VectorMVector, 
            f_tmp::VectorMVector, f::VectorMVector, y1, y2, error, 
            J::MatrixMMatrix, linear_cache, finitediff_cache,
            jacobian_config, dy_dt_wrap!, args...) where {T <: AbstractFloat, F}

    fixed_runge_kutta_step!(method, iteration, y, t[1], dt[1], dy_dt!, dy, y_tmp, f_tmp,
                            J, linear_cache, finitediff_cache, jacobian_config,
                            dy_dt_wrap!
                            )
    y .= y_tmp
    nothing
end
