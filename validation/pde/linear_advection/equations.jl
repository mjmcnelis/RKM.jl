
function dy_dt_central!(f, y, t; p, kwargs...)
    a, dx, dt = p
    N = length(y)
    for i in 1:N
        m = max(i-1, 1) # BC: y[0] = y[1]
        p = min(i+1, N) # BC: y[N+1] = y[N]
        ym, yp = y[m], y[p]
        Fp = a*yp
        Fm = a*ym
        f[i] = -(Fp - Fm)/(2.0*dx)
    end
    nothing
end

function dy_dt_lax_friedichs!(f, y, t; p, kwargs...)
    a, dx, _ = p
    N = length(y)
    for i in 1:N
        m = max(i-1, 1) # BC: y[0] = y[1]
        p = min(i+1, N) # BC: y[N+1] = y[N]
        ym, yc, yp = y[m], y[i], y[p]
        Fp = a*yp
        Fm = a*ym
        f[i] = -(Fp - Fm)/(2.0*dx) + (yp - 2.0*yc + ym)/(2.0*dt)
    end
    nothing
end

function dy_dt_upwind!(f, y, t; p, kwargs...)
    a, dx, _ = p
    N = length(y)
    for i in 1:N
        m = max(i-1, 1) # BC: y[0] = y[1]
        p = min(i+1, N) # BC: y[N+1] = y[N]
        ym, yc, yp = y[m], y[i], y[p]
        Fp = a > 0.0 ? a*yc : a*yp
        Fm = a > 0.0 ? a*ym : a*yc
        f[i] = -(Fp - Fm)/dx
    end
    nothing
end

function dy_dt_rusanov!(f, y, t; p, kwargs...)
    a, dx, _ = p
    N = length(y)
    for i in 1:N
        m = max(i-1, 1) # BC: y[0] = y[1]
        p = min(i+1, N) # BC: y[N+1] = y[N]
        ym, yc, yp = y[m], y[i], y[p]
        Fp = a*yp
        Fm = a*ym
        f[i] = -(Fp - Fm)/(2.0*dx) + (yp - 2.0*yc + ym)*a/(2.0*dx)
    end
    nothing
end

function gauss(x)
    exp(-(x - 3)^2)
end
