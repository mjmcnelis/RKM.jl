
function dy_dt!(f, y, t; p, kwargs...)
    a, dx = p
    N = length(y)
    for i in 1:N
        m = max(i-1, 1) # BC: y[0] = y[1]
        p = min(i+1, N) # BC: y[N+1] = y[N]
        ym, yc, yp = y[m], y[i], y[p]
        f[i] = (yp - 2.0*yc + ym)*a/dx^2
    end
    nothing
end

# TODO: rework exact solution under finite x-interval and NBC
function gauss(x, t; p, t0)
    a = p[1]
    return sqrt(t0/t) * exp(-x^2/(4.0*a*t))
end

# OrdinaryDiffEq
function diffeq!(f, y, p, t)
    a, dx = p
    N = length(y)
    for i in 1:N
        m = max(i-1, 1) # BC: y[0] = y[1]
        p = min(i+1, N) # BC: y[N+1] = y[N]
        ym, yc, yp = y[m], y[i], y[p]
        f[i] = (yp - 2.0*yc + ym)*a/dx^2
    end
    nothing
end