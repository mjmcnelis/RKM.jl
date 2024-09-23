
function dy_dt!(f, y; p, kwargs...)
    for i in eachindex(y)
        f[i] = (y[i] + p[i]) * (1.0 - p[i] - y[i])
    end
    nothing
end

function y_exact(t; N = 1)
    y_ex = BigFloat[]
    t = BigFloat(t)
    for i = 1:N
        p = 0.5 - 0.25*(i-1.0)/(N-1.0+eps(1.0))
        push!(y_ex, exp(t) / (1.0 + exp(t)) - p)
    end
    y_ex
end

# for OrdinaryDiffEq (StaticArray version)
function f_ord_static(y, p, t)
    SA[((y[i] + p[i])*(1.0 - p[i] - y[i]) for i in eachindex(y))...]
end

# for OrdinaryDiffEq (normal version)
function f_ord(f, y, p, t)
    for i in eachindex(y)
        f[i] = (y[i] + p[i]) * (1.0 - p[i] - y[i])
    end
    nothing
end
