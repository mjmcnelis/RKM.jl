
# for RKM 
function dy_dt!(f, t, y)
    f[1] = - t * (y[1] - 1.0)
    return nothing 
end

# for OrdinaryDiffEq
function f_ord(f, y, p, t)
    f[1] = - t * (y[1] - 1.0)
    return nothing 
end

# exact solution 
function y_exact(t)
    t = BigFloat(t)
    return BigFloat[1.0 + exp(-t^2/2.0)]
end