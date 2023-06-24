# for RKM
function dy_dt!(f, y; kwargs...)
    f[1] = y[1]^2*(1.0 - y[1])
    return nothing
end

# for OrdinaryDiffEq
function dy_dt!(f, y, p, t)
    f[1] = y[1]^2*(1.0 - y[1])
    return nothing
end
