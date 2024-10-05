
# for RKM
function dy_dt!(f, y, t; p, kwargs...)
    ω = p[1]
    f[1] = y[2]
    f[2] = -ω^2*y[1]
    return nothing
end

# for OrdinaryDiffEq
function dy_dt!(f, y, p, t)
    ω = p[1]
    f[1] = y[2]
    f[2] = -ω^2*y[1]
    return nothing
end
