# for RKM
function dy_dt!(f, y; p, kwargs...)
    σ, ρ, β = p
    f[1] = σ*(y[2] - y[1])
    f[2] = y[1]*(ρ - y[3]) - y[2]
    f[3] = y[1]*y[2] - β*y[3]
    return nothing
end

# for OrdinaryDiffEq
function dy_dt!(f, y, p, t)
    σ, ρ, β = p
    f[1] = σ*(y[2] - y[1])
    f[2] = y[1]*(ρ - y[3]) - y[2]
    f[3] = y[1]*y[2] - β*y[3]
    return nothing
end
