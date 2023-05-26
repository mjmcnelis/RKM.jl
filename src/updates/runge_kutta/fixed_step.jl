
function evolve_one_time_step!(method::RungeKutta, iteration::Explicit,
             adaptive::Fixed, ::Controller, FE::MVector{1,Int64}, 
             y::VectorMVector, t::VectorMVector{1,T}, dt::VectorMVector{2,T}, 
             dy_dt!::F, dy::MatrixMMatrix, y_tmp::VectorMVector, 
             f_tmp::VectorMVector, args...) where {T <: AbstractFloat, F}
          
    # note: I can comment this out and loop i = 1:stages w/o allocating
    dy_dt!(f_tmp, t[1], y)                              # evaluate first stage at (t,y)
    @.. dy[:,1] = dt[1] * f_tmp

    runge_kutta_step!(method, iteration, y, t[1], dt[1], dy_dt!, dy, y_tmp, f_tmp)
    @.. y = y_tmp                                       # get iteration

    add_function_evaluations!(FE, iteration, adaptive, method)
    return nothing
end

function evolve_one_time_step!(method::RungeKutta, iteration::DiagonalImplicit,
             adaptive::Fixed, ::Controller, FE::MVector{1,Int64},
             y::VectorMVector, t::VectorMVector{1,T}, dt::VectorMVector{2,T},
             dy_dt!::F, dy::MatrixMMatrix, y_tmp::VectorMVector, 
             f_tmp::VectorMVector, f::VectorMVector, y1, y2, error, 
             J::MatrixMMatrix, linear_cache, dy_dt_wrap!,
             stage_finder::ImplicitStageFinder) where {T <: AbstractFloat, F}

    runge_kutta_step!(method, iteration, y, t[1], dt[1], dy_dt!, dy, y_tmp, 
                      f_tmp, J, linear_cache, dy_dt_wrap!, stage_finder)
    @.. y = y_tmp
    return nothing
end