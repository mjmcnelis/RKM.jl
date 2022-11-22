
function dy_dt!(f, t, y)
    for i in eachindex(y)
        a = 0.5 - (i - 1.0)/4.0
        f[i] = (y[i] + a) * (1.0 - a - y[i])
    end
    nothing
end

function jacobian!(J, t, y)
    for i in eachindex(y)
        a = 0.5 - (i - 1.0)/4.0
        J[i,i] = 1.0 - 2.0*(y[i] + a)
    end
    nothing 
end
