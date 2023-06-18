
get_a(i, N) = N == 1 ? 0.5 : 0.5 - 0.25*(i-1.0)/(N-1.0)

function dy_dt!(f, y; kwargs...)
    N = length(y)
    for i in eachindex(y)
        a = get_a(i, N)
        f[i] = (y[i] + a) * (1.0 - a - y[i])
    end
    nothing
end

# pass model_parameters = [0.5, 0.25]
# function dy_dt!(f, y; p, kwargs...)
#     N = length(y)
#     for i in eachindex(y)
#         a = p[i]
#         f[i] = (y[i] + a) * (1.0 - a - y[i])
#     end
#     nothing
# end

function y_exact(t; N = 1)
    y_ex = BigFloat[]
    t = BigFloat(t)
    for i = 1:N
        a = N == 1 ? 0.5 : 0.5 - 0.25*(i-1.0)/(N-1.0)
        push!(y_ex, exp(t) / (1.0 + exp(t)) - a)
    end
    y_ex
end

# for OrdinaryDiffEq (StaticArray version)
function f_ord_static(y, p, t)
    N = length(y)
    SA[((y[i] + get_a(i,N))*(1.0 - get_a(i,N) - y[i]) for i in eachindex(y))...]
end

# for OrdinaryDiffEq (normal version)
function f_ord(f, y, p, t)
    N = length(y)
    for i in eachindex(y)
        a = get_a(i, N)
        f[i] = (y[i] + a) * (1.0 - a - y[i])
    end
    nothing
end
