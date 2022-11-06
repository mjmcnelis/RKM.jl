
function dy_dt!(f, t, y)
    a = 0.5
    b = 0.25
    f[1] = (y[1] + a) * (1.0 - a - y[1])
    f[2] = (y[2] + b) * (1.0 - b - y[2])
    nothing
end
