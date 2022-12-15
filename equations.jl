
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

# for OrdinaryDiffEq (StaticArray version)
function fp(y, p, t)
    SA[(y[1] + 0.5) * (0.5 - y[1]),
       (y[2] + 0.25) * (0.75 - y[2])]
end
# for OrdinaryDiffEq (normal version)
# function fp(f, y, p, t)
#     f[1] = (y[1] + 0.5) * (0.5 - y[1])
#     f[2] = (y[2] + 0.25) * (0.75 - y[2])
#     nothing
# end