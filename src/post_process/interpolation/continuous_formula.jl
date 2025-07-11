
@muladd function interpolate_solution(interpolator::ContinuousFormula, sol::Solution,
                                      method::ODEMethod, precision::Type{T};
                                      dt_dense::Float64) where T <: AbstractFloat

    t = sol.t
    y = sol.y
    dy = sol.dy
    dimensions = sol.dimensions

    if isempty(y)
        error("Original solution set is empty, set save_solution = true")
    end

    ω = method.ω
    stages = method.stages
    reconstructor = method.reconstructor

    if isempty(ω)
        error("$reconstructor has no coefficients for ContinuousFormula \
               (use CubicHermite instead)")
    end

    # order of continuous formula
    q = size(ω, 2)
    q_vect = 1:q |> SVector{q, Int64}

    @info "Generating order-$q continuous output for $reconstructor"

    # get dimensions
    nt = length(t)
    ny = dimensions[1]

    # initialize intermediate cache
    y_interp = zeros(precision, ny)
    y_interp .= view(y, 1:ny)

    # dense time
    t0 = t[1]
    tf = t[end]
    nt_dense = 1 + floor(Int64, Float64((tf - t0)/dt_dense))
    t_dense = range(t0, tf, nt_dense)

    # dense state vector
    y_dense = Vector{precision}()
    sizehint!(y_dense, ny*nt_dense)
    append!(y_dense, y_interp)

    # loop over original solution
    for n in 1:nt-1
        idxs_y = (1 + (n-1)*ny):n*ny
        idxs_dy = (1 + (n-1)*ny*stages):n*ny*stages

        y1 = view(y, idxs_y)
        dy1 = reshape(view(dy, idxs_dy), ny, stages)

        t1 = t[n]
        t2 = t[n+1]
        Δt = t2 - t1

        nL = 2 + floor(Int64, Float64((t1 - t0)/dt_dense))
        nR = 1 + floor(Int64, Float64((t2 - t0)/dt_dense))

        # loop over dense points
        for m in nL:min(nR, nt_dense)
            t_interp = t_dense[m]
            θ = (t_interp - t1) / Δt
            θ_pow = θ.^q_vect           # [θ, θ^2, ..., θ^q]
            bθ = ω * θ_pow              # bθ_j = ω_jk * θ^k

            # continuous output formula
            @.. y_interp = y1
            for j in 1:stages
                dy_stage = view(dy1,:,j)
                @.. y_interp = y_interp + bθ[j]*dy_stage
            end

            # store interpolated solution
            append!(y_dense, y_interp)
        end
    end
    # reshape dense solution
    y_dense = reshape(y_dense, ny, nt_dense) |> transpose

    return t_dense, y_dense
end