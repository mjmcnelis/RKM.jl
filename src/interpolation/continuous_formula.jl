
@muladd function interpolate_solution(interpolator::ContinuousFormula, sol::Solution,
                                      precision::Type{T} = Float64;
                                      dt_dense::Float64) where T <: AbstractFloat

    @unpack t, y, f, dimensions = sol

    if isempty(y)
        @warn "Original solution is empty, dense output will also be empty..."
        return get_solution(sol)
    end

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

    # TMP:
    #--------------------------
    # get continuous coefficients ω from method
    ω = SMatrix{3, 4, precision, 12}(
        1, -3//2, 2//3,
        0, 1, -2//3,
        0, 1, -2//3,
        0, -1//2, 2//3
    ) |> transpose

    s = size(ω, 1)  # number of stages
    q = size(ω, 2)  # order of continuous formula
    #--------------------------

    bθ = MVector{s, precision}(zeros(s))
    θ_pow = MVector{p, precision}(zeros(q))
    q_vect = SVector{q, Int64}(1:q)

    # loop over original solution
    for n in 1:nt-1
        idxs_y = (1 + (n-1)*ny):n*ny

        y1 = view(y, idxs_y)

        # TODO: get intermediate stages instead
        # f1 = view(f, idxs_1)
        # f2 = view(f, idxs_2)

        t1 = t[n]
        t2 = t[n+1]
        Δt = t2 - t1

        nL = 2 + floor(Int64, Float64((t1 - t0)/dt_dense))
        nR = 1 + floor(Int64, Float64((t2 - t0)/dt_dense))

        # loop over dense points
        for m in nL:min(nR, nt_dense)
            t_interp = t_dense[m]
            θ = (t_interp - t1) / Δt
            θ_pow .= θ.^q_vect          # [θ, θ^2, ..., θ^q]
            bθ .= ω * θ_pow             # bθ_j = ω_jk * θ^k

            # continuous output formula
            @.. y_interp = y1
            for j in 1:s
                dy_stage = view(dy,:,j)
                @.. y_interp = y_interp + bθ[j]*dy_stage
            end

            # store interpolated solution
            append!(y_dense, y_interp)
        end
    end
    # reshape dense solution
    # y_dense = reshape(y_dense, ny, nt_dense) |> transpose

    return t_dense, y_dense
end