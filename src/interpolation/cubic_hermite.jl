
function interpolate_solution(interpolator::CubicHermite, sol::Solution,
                              method::ODEMethod, precision::Type{T};
                              dt_dense::Float64) where T <: AbstractFloat

    @unpack t, y, f, dimensions = sol

    if isempty(y)
        error("Original solution set is empty, set save_solution = true")
    end

    @info "Generating cubic Hermite interpolation"

    # get dimensions
    nt = length(t)
    ny = dimensions[1]

    # initialize intermediate cache
    y_interp = zeros(precision, ny)
    y_interp .= view(y, 1:ny)

    # dense time
    # TODO: pass t_dense directly (but need to fix nL/nR)
    t0 = t[1]
    tf = t[end]
    nt_dense = 1 + floor(Int64, Float64((tf - t0)/dt_dense))
    t_dense = range(t0, tf, nt_dense)

    # dense state vector
    y_dense = Vector{precision}()
    sizehint!(y_dense, ny*nt_dense)
    append!(y_dense, y_interp)

    # TODO: can I get dense derivatives easily?

    # loop over original solution
    for n in 1:nt-1
        idxs_1 = (1 + (n-1)*ny):n*ny
        idxs_2 = (1 + n*ny):(n+1)*ny

        y1 = view(y, idxs_1)
        y2 = view(y, idxs_2)

        f1 = view(f, idxs_1)
        f2 = view(f, idxs_2)

        t1 = t[n]
        t2 = t[n+1]
        Δt = t2 - t1

        nL = 2 + floor(Int64, Float64((t1 - t0)/dt_dense))
        nR = 1 + floor(Int64, Float64((t2 - t0)/dt_dense))

        # loop over dense points
        for m in nL:min(nR, nt_dense)
            t_interp = t_dense[m]
            θ = (t_interp - t1) / Δt
            θ2 = θ * θ
            θ3 = θ2 * θ

            # Hermite interpolation formula
            @.. y_interp = (2.0*(y1-y2) + (f1+f2)*Δt)*θ3 +
                           (3.0*(y2-y1) - (2.0*f1+f2)*Δt)*θ2 + f1*Δt*θ + y1

            # store interpolated solution
            append!(y_dense, y_interp)
        end
    end
    # reshape dense solution
    y_dense = reshape(y_dense, ny, nt_dense) |> transpose

    return t_dense, y_dense
end
