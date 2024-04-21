
abstract type Interpolator end

struct NoInterpolator <: Interpolator end

@kwdef struct HermiteInterpolator <: Interpolator
    dt_save::Float64
    # t_range::StepRangeLen
end

#=
function HermiteInterpolator(; dt_save::Float64)
    # t0 = 0.0
    # tf = 1.0

    # nt = 1 + floor(Int, (tf - t0)/dt_save)
    # t_range = range(t0, t0 + (nt-1)*dt_save, nt)

    # t_range = range()

    # @show t_range
    # q()
    # include save_counter so can index t_range correctly

    return HermiteInterpolator(dt_save, t_range)
end
=#

function interpolate_solution!(::NoInterpolator, sol::Solution,
                               update_cache::RKMCache, t::Vector{T},
                               t0::T, tf::T) where T <: AbstractFloat
    @unpack y = update_cache

    append!(sol.y, y)
    append!(sol.t, t[1])
    return nothing
end

function interpolate_solution!(interp::HermiteInterpolator, sol::Solution,
                               update_cache::RKMCache, t::Vector{T},
                               t0::T, tf::T) where T <: AbstractFloat
    @unpack dt_save = interp
    @unpack y, y1, y2, y_tmp, f, f_tmp = update_cache

    # this is the first thing I should take care of
    # store previous solution in y_tmp, what about f_tmp?
    # q()

    nt = 1 + floor(Int, (tf - t0)/dt_save)
    t_range = range(t0, t0 + (nt-1)*dt_save, nt)

    # @show t_range
    t_tmp = t_range[1]
    # @show t_tmp
    # q()
    # q()

    # should call interpolate after one time step
    # the number of times I interpolate depends on how many saved points I pass

    # @show y1
    # q()

    return nothing
end