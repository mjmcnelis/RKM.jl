
function RKM.norm(v::Array{DoubleFloat{T}, N}, p::Real = 2.0) where {N, T<:IEEEFloat}
    isempty(v) && return zero(DoubleFloat{T})

    if isinf(p)
        return signbit(p) ? minimum(abs, v) : maximum(abs, v)
    elseif p == 2
        return sqrt(sum(abs2, v))
    else
        # note: had more compact version (but deleted file...)
        vp = DoubleFloat{T}(0.0)
        for i in eachindex(v)
            vp += abs(v[i])^p
        end
        r = inv(DoubleFloat{T}(p))
        return vp^r
    end
end
