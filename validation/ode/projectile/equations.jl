
function dy_dt!(f, y, t; p, kwargs...)
    γ, g = p
    f[1] = y[2]
    f[2] = -g - γ*y[2]
    return nothing
end

function y_exact(t; y0, t0, p)
    h0, v0 = y0 # initial height, velocity
    γ, g = p    # damping coefficient and gravitational acceleration

    y = zeros(BigFloat, 2)
    t = BigFloat(t)

    if γ == 0.0
        y[1] = h0 + v0*(t-t0) - 0.5*g*(t-t0)^2  # height
        y[2] = v0 - g*(t-t0)                    # velocity
    else
        y[1] = h0 - g*(t-t0)/γ + (v0 + g/γ)/γ*(1.0 - exp(-γ*(t-t0)))
        y[2] = -g/γ + (v0 + g/γ)*exp(-γ*(t-t0))
    end
    return y
end