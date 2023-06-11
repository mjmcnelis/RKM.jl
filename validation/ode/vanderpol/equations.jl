# for RKM
function dy_dt!(f, t, y)
    mu = 1000.0
    f[1] = y[2] 
    f[2] = mu*(1.0 - y[1]^2)*y[2] - y[1]
    return nothing 
end

# for OrdinaryDiffEq 
function dy_dt!(f, y, p, t)
    mu = 1000.0
    f[1] = y[2] 
    f[2] = mu*(1.0 - y[1]^2)*y[2] - y[1]
    return nothing 
end
