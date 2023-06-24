
using Revise
using RKM
using LinearAlgebra
using Plots
using UnPack
plotly()

# initial condition
x = LinRange(-5, 5, 101) |> collect
const a = 0.25            # diffusion constant
const dx = x[2] - x[1]
const dt = 0.01           # just adjust this for courant number
const C = 2*a*dt/dx^2
const t0 = 1.0
N = 300

# TODO: rework exact solution under finite x-interval and NBC
function gauss(x, t)
    sqrt(t0/t) * exp(-x^2/(4*a*t))
end

function dy_dt!(f, t, y)
    L = length(y)
    for i in 1:L
        m = max(i-1, 1) # BC: y[0] = y[1]
        p = min(i+1, L) # BC: y[L+1] = y[L]
        ym, yc, yp = y[m], y[i], y[p]

        f[i] = (yp - 2.0*yc + ym)*a/dx^2
    end
    nothing
end

adaptive   = Fixed()
method     = Euler1()
t_range    = TimeRange(; t0, tf = 10.0, dt0 = dt)
parameters = Parameters(; adaptive, method, t_range)

@unpack t0, dt0 = t_range
y0 = gauss.(x, t0)
@show C

@time sol = evolve_ode(y0, dy_dt!; parameters)

plt = plot(x, y0, label = "t = $t0", color = "indianred", linewidth = 2,
           size = (900, 600), ylims = (-0.5, 1.3),
           ylabel = "u", yguidefontsize = 14, ytickfontsize = 12,
           xlabel = "t", xguidefontsize = 14, xtickfontsize = 12,
           legend = :outertopright, legendfontsize = 12)
for i = 1:3
    t = t0 + N*i*dt0
    y = sol.y[1 + N*i]
    plot!(x, y, color = "indianred", linewidth = 2, label = "t = $t")
end

plot!(x, y0, label = "t = $t0 (exact)", color = "black", linewidth = 1, line = :dash)
for i = 1:3
    t = t0 + N*i*dt0
    y_exact = gauss.(x, t)
    plot!(x, y_exact, color = "black", linewidth = 1, line = :dash, label = "t = $t (exact)")
end
display(plt)

println("\ndone")

