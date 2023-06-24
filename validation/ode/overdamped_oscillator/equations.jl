# for RKM
function dy_dt!(f, y; p, kwargs...)
    γ, ω = p
    f[1] = y[2]
    f[2] = -γ*y[2] - ω^2*y[1]
    return nothing
end

# for OrdinaryDiffEq
function dy_dt!(f, y, p, t)
    γ, ω = p
    f[1] = y[2]
    f[2] = -γ*y[2] - ω^2*y[1]
    return nothing
end
