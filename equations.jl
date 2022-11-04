
function dy_dt!(f, t, y)
    A = 0.5
    f[1] = (y[1] + A) * (1.0 - A - y[1])
    # f .= (y .+ A) .* (1.0 .- A .- y)
    nothing
end
