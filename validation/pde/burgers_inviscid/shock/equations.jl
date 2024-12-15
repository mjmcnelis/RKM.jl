# flux term in Burgers inviscid equation
F(y) = y^2/2.0

function dy_dt_central!(f, y, t; p, kwargs...)
    dx, _ = p
    N = length(y)
    for i in 1:N
        m = max(i-1, 1)                 # NBC: y[0] = y[1]
        p = min(i+1, N)                 # NBC: y[N+1] = y[N]
        ym, yp = y[m], y[p]

        Fp = F(yp)                      # numerical fluxes
        Fm = F(ym)

        f[i] = -(Fp - Fm)/(2.0*dx)
    end
    nothing
end

function dy_dt_rusanov!(f, y, t; p, kwargs...)
    dx, _ = p
    N = length(y)
    for i in 1:N
        m = max(i-1, 1)                 # NBC: y[0] = y[1]
        p = min(i+1, N)                 # NBC: y[N+1] = y[N]
        ym, yc, yp = y[m], y[i], y[p]

        aR = max(abs(yp), abs(yc))      # local propagation speeds
        aL = max(abs(yc), abs(ym))

        Fp = F(yp) - aR*(yp - yc)       # numerical fluxes
        Fm = F(ym) - aL*(yc - ym)

        f[i] = -(Fp - Fm)/(2.0*dx)
    end
    nothing
end

function shock(x)
    x < 0.0 ? 1.0 :
    x > 0.0 ? 0.0 : 0.5
end
