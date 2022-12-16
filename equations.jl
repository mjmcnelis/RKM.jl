
function dy_dt!(f, t, y)
    N = length(y)
    for i in eachindex(y)
        a = N == 1 ? 0.5 : 0.5 - 0.25*(i-1.0)/(N-1.0)
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
get_a(i, N) = N == 1 ? 0.5 : 0.5 - 0.25*(i-1.0)/(N-1.0)

function fp(y, p, t)
    N = length(y)
    SA[((y[i] + get_a(i,N))*(1.0 - get_a(i,N) - y[i]) for i = 1:N)...]
end
# for OrdinaryDiffEq (normal version)
# function fp(f, y, p, t)
#     f[1] = (y[1] + 0.5) * (0.5 - y[1])
#     f[2] = (y[2] + 0.25) * (0.75 - y[2])
#     nothing
# end