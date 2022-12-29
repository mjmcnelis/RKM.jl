
function dy_dt!(f, t, y)
    N = length(y)
    for i in 1:N
        a = N == 1 ? 0.5 : 0.5 - 0.25*(i-1.0)/(N-1.0)
        f[i] = (y[i] + a) * (1.0 - a - y[i])
    end
    nothing
end

function y_exact(t, N)
    y_ex = Float64[]
    for i = 1:N
        a = N == 1 ? 0.5 : 0.5 - 0.25*(i-1.0)/(N-1.0)
        push!(y_ex, exp(t) / (1.0 + exp(t)) - a)
    end
    y_ex
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
#     N = length(y)
#     for i in 1:N
#         a = N == 1 ? 0.5 : 0.5 - 0.25*(i-1.0)/(N-1.0)
#         f[i] = (y[i] + a) * (1.0 - a - y[i])
#     end
#     nothing
# end
