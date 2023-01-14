
# for RKM 
function dy_dt!(f, t, y)
    f[1] = y[2]
    f[2] = -ω^2*y[1]
    return nothing 
end

# for OrdinaryDiffEq 
function fp(f, y, p, t)
    f[1] = y[2]
    f[2] = -ω^2*y[1]
    return nothing 
end
