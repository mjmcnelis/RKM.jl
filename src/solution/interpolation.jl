
abstract type Interpolator end
abstract type DenseInterpolator <: Interpolator end

struct NoInterpolator <: Interpolator end

struct HermiteInterpolator{SRL} <: DenseInterpolator where SRL <: StepRangeLen
    dt_save::Float64
    t_range::SRL
    nt::Int64
    t0::Float64
end

function HermiteInterpolator(; dt_save::Float64)
    # note: nt, t0 and t_range are dummy values
    nt = 2
    t0 = 0.0
    t_range = range(t0, 1.0, nt)
    return HermiteInterpolator(dt_save, t_range, nt, t0)
end

reconstruct_interpolator(interpolator::NoInterpolator, args...) = interpolator

function reconstruct_interpolator(interpolator::HermiteInterpolator,
                                  t0::T, tf::T) where T <: AbstractFloat
    @unpack dt_save = interpolator
    nt = 1 + floor(Int64, Float64((tf - t0)/dt_save))
    t_range = range(t0, tf, nt)

    return HermiteInterpolator(dt_save, t_range, nt, Float64(t0))
end

function interpolate_solution!(::NoInterpolator, sol::Solution, update_cache::RKMCache,
                               t::Vector{T}) where T <: AbstractFloat
    @unpack y_tmp, S = update_cache

    append!(sol.y, y_tmp)
    append!(sol.t, t[1])

    # TODO: skip if not getting sensitivities
    append!(sol.S, S)       # use S_tmp instead once make it
    return nothing
end

function interpolate_solution!(interpolator::HermiteInterpolator,
                               sol::Solution, update_cache::RKMCache,
                               t::Vector{T}) where T <: AbstractFloat

    @unpack dt_save, t_range, t0, nt = interpolator

    t_curr, t_prev = t
    Δt = t_curr - t_prev

    nL = 2 + floor(Int64, Float64((t_prev - t0)/dt_save))
    nR = 1 + floor(Int64, Float64((t_curr - t0)/dt_save))

    @unpack dy, y, y_tmp, f, f_tmp = update_cache
    y1, y2, f1, f2 = y, y_tmp, f, f_tmp
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
