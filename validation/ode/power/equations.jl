# for RKM
function dy_dt!(f, y, t; p, kwargs...)
    n = p[1]
    f[1] = n * y[1]^((n-1.0)/n)
    return nothing
end

# exact solution
function y_exact(t; p)
    t = BigFloat(t)
    n = p[1]
    return t^n
end