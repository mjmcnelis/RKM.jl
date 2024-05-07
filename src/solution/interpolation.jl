
abstract type Interpolator end
abstract type DenseInterpolator <: Interpolator end

struct NoInterpolator <: Interpolator end

struct HermiteInterpolator{SRL} <: DenseInterpolator where SRL <: StepRangeLen
    dt_save::Float64
    t_range::SRL
    nt::Int64
end

function HermiteInterpolator(; dt_save::Float64)
    # note: nt and t_range are dummy values
    nt = 2
    t_range = range(0, 1, nt)
    return HermiteInterpolator(dt_save, t_range, nt)
end

reconstruct_interpolator(interpolator::NoInterpolator, args...) = interpolator

function reconstruct_interpolator(interpolator::HermiteInterpolator,
                                  t0::T, tf::T) where T <: AbstractFloat
    @unpack dt_save = interpolator
    nt = 1 + floor(Int64, Float64((tf - t0)/dt_save))
    t_range = range(t0, tf, nt)

    return HermiteInterpolator(dt_save, t_range, nt)
end

function interpolate_solution!(::NoInterpolator, sol::Solution, update_cache::RKMCache,
                               t::Vector{T}, args...) where T <: AbstractFloat
    @unpack y = update_cache

    append!(sol.y, y)
    append!(sol.t, t[1])
    return nothing
end

function interpolate_solution!(interpolator::HermiteInterpolator, sol::Solution,
                               update_cache::RKMCache, t::Vector{T},
                               t0::T) where T <: AbstractFloat

    # should call interpolate after one time step taken
    # the number of times I interpolate depends on how many saved points I pass

    # for FSAL methods, I already have f = y_n and f_tmp = y_{n+1}
    # but in general I would have to do a function evaluation for f_tmp
    # doesn't mean it'll cost me an extra FE, I just have to save previous evaluation
    # ode_wrap!(f_tmp, t[1], y)
    # in next time step I set f .= f_tmp

    @unpack dt_save, t_range, nt = interpolator

    t_curr, t_prev = t
    Δt = t_curr - t_prev

    nL = 2 + floor(Int64, Float64((t_prev - t0)/dt_save))
    nR = 1 + floor(Int64, Float64((t_curr - t0)/dt_save))

    @unpack dy, y, y_tmp, f, f_tmp = update_cache
    y1, y2, f1, f2 = y_tmp, y, f, f_tmp
    y_interp = view(dy,:,1)

    for n in nL:min(nR, nt)
        t_interp = t_range[n]
        Θ = (t_interp - t_prev) / Δt
        Θ2 = Θ * Θ
        Θ3 = Θ2 * Θ

        @.. y_interp = (2*(y1-y2) + (f1+f2)*Δt)*Θ3 +
                       (3*(y2-y1) - (2*f1+f2)*Δt)*Θ2 + f1*Δt*Θ + y1

        append!(sol.y, y_interp)
        append!(sol.t, t_interp)
    end

    return nothing
end
