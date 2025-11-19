function dy_dt!(f, y, t; p, abstract_params)
    # note: here dx is not a sensitivity parameter
    a = p[1]
    dx = abstract_params
    N = length(y)
    for i in 1:N
        # Dirichlet boundary conditions: y[0] = 0.0, y[N+1] = 0.0
        ym = i == 1 ? 0.0 : y[i-1]
        yp = i == N ? 0.0 : y[i+1]
        yc = y[i]
        f[i] = (yp - 2.0*yc + ym)*a/dx^2
    end
    nothing
end

function gauss(x, t; p, t0)
    a = p[1]
    return sqrt(t0/t) * exp(-x^2/(4.0*a*t))
end