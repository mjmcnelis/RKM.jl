
abstract type Interpolator end
abstract type DenseInterpolator <: Interpolator end

struct NoInterpolator <: Interpolator end

@kwdef struct HermiteInterpolator <: DenseInterpolator
    dt_save::Float64
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

    # should call interpolate after one time step taken
    # the number of times I interpolate depends on how many saved points I pass

    # for FSAL methods, I already have f = y_n and f_tmp = y_{n+1}
    # but in general I would have to do a function evaluation for f_tmp
    # doesn't mean it'll cost me an extra FE, I just have to save previous evaluation
    # ode_wrap!(f_tmp, t[1], y)
    # in next time step I set f .= f_tmp

    dt_save = T(interp.dt_save)

    t_curr, t_prev = t
    Δt = t_curr - t_prev

    nL = 2 + floor(Int64, Float64((t_prev - t0)/dt_save))
    nR = 1 + floor(Int64, Float64((t_curr - t0)/dt_save))

    # shouldn't need to compute nt, t_range every time
    nt = 1 + floor(Int64, Float64((tf - t0)/dt_save))
    # t_range = range(t0, t0 + (nt-1)*dt_save, nt)
    t_range = range(t0, tf, nt)

    # @show nL, nR, nt t_prev, t_curr
    # q()

    @unpack dy, y, y_tmp, f, f_tmp = update_cache
    y1, y2, f1, f2 = y_tmp, y, f, f_tmp
    y_interp = view(dy,:,1)

    for n in nL:min(nR, nt)
        Θ = (t_range[n] - t_prev) / Δt
        Θ2 = Θ * Θ
        Θ3 = Θ2 * Θ

        @.. y_interp = (2*(y1-y2) + (f1+f2)*Δt)*Θ3 +
                       (3*(y2-y1) - (2*f1+f2)*Δt)*Θ2 + f1*Δt*Θ + y1

        append!(sol.y, y_interp)
        append!(sol.t, t_range[n])
    end

    return nothing
end
