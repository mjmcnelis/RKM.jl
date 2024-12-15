function dy_dt!(f, y, t; p, kwargs...)
    a, dx, dt = p
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