
function evolve_one_time_step!(method::RungeKutta, iteration::Iteration,
             adaptive::Fixed, controller::Controller, FE::MVector{1,Int64},
             y::VectorMVector, t::VectorMVector{1,T}, dt::VectorMVector{2,T},
             dy_dt!::F, dy::MatrixMMatrix, y_tmp::VectorMVector, 
             f_tmp::VectorMVector, f::VectorMVector, y1, y2, error, 
             J::MatrixMMatrix, linear_cache, dy_dt_wrap!::ODEWrapper,
             stage_finder::ImplicitStageFinder) where {T <: AbstractFloat, F}
    # note: for explicit, can comment this and loop i = 1:stages w/o allocating
    if iteration isa Explicit
        dy_dt!(f_tmp, t[1], y)                          # evaluate first stage at (t,y)
        @.. dy[:,1] = dt[1] * f_tmp
        FE[1] += 1
    end

    runge_kutta_step!(method, iteration, y, t[1], dt[1], dy_dt!, dy, y_tmp, 
                      f_tmp, FE, J, linear_cache, dy_dt_wrap!, stage_finder)
    @.. y = y_tmp
    return nothing
end