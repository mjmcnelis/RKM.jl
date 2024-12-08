# for RKM
function dy_dt!(f, y, t; p, kwargs...)
    f[1] = p[1] * y[1]
    return nothing
end

# exact solution
function y_exact(t; p)
    t = BigFloat(t)
    return exp(p[1]*t)
end