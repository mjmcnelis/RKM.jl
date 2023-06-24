
# for RKM
function dy_dt!(f, y; kwargs...)
    k1 = 0.04
    k2 = 3.0e7
    k3 = 1.0e4
    f[1] = -k1*y[1] + k3*y[2]*y[3]
    f[2] = k1*y[1] - k2*y[2]^2 - k3*y[2]*y[3]
    f[3] = k2*y[2]^2
    return nothing
end

# for OrdinaryDiffEq
function dy_dt!(f, y, p, t)
    k1 = 0.04
    k2 = 3.0e7
    k3 = 1.0e4
    f[1] = -k1*y[1] + k3*y[2]*y[3]
    f[2] = k1*y[1] - k2*y[2]^2 - k3*y[2]*y[3]
    f[3] = k2*y[2]^2
    return nothing
end
