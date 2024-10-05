# for RKM
function dy_dt!(f, y, t; p, kwargs...)
    α, β, γ, δ = p
    f[1] = (α - β*y[2])*y[1]
    f[2] = (δ*y[1] - γ)*y[2]
    return nothing
end

# for OrdinaryDiffEq
function dy_dt!(f, y, p, t)
    α, β, γ, δ = p
    f[1] = (α - β*y[2])*y[1]
    f[2] = (δ*y[1] - γ)*y[2]
    return nothing
end
