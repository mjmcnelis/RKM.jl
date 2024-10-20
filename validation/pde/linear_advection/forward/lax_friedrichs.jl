
using Revise
using RKM
using LinearAlgebra
using Plots
using UnPack
plotly()

function gauss(x)
    exp(-(x - 3)^2)
end

# initial condition
x = LinRange(0, 12, 121) |> collect
const a = 1.0
const dx = x[2] - x[1]
const dt = 0.05           # just adjust this for courant number
const C = a*dt/dx
N = 40

F(y) = a*y

function dy_dt!(f, t, y)
    L = length(y)
    for i in 1:L
        m = max(i-1, 1) # BC: y[0] = y[1]
        p = min(i+1, L) # BC: y[L+1] = y[L]
        ym, yc, yp = y[m], y[i], y[p]

        f[i] = -(F(yp) - F(ym))/(2.0*dx) + (yp - 2.0*yc + ym)/(2.0*dt)
    end
    nothing
end

adaptive   = Fixed()
method     = Euler1()
t_range    = TimeRange(; t0 = 0.0, tf = 6.0, dt0 = dt)
options = SolverOptions(; adaptive, method, t_range)

@unpack t0, dt0 = t_range
y0 = gauss.(x)
@show C

@time sol = evolve_ode(y0, dy_dt!; parameters)
t, y = get_solution(sol)

plt = plot(x, y0, label = "t = 0", color = "indianred", linewidth = 2,
           size = (900, 600), ylims = (-0.5, 1.3),
           ylabel = "u", yguidefontsize = 14, ytickfontsize = 12,
           xlabel = "x", xguidefontsize = 14, xtickfontsize = 12,
           legend = :outertopright, legendfontsize = 12, dpi = 200)
for i = 1:3
    t = t0 + N*i*dt0
    plot!(x, y[1+N*i, :], color = "indianred", linewidth = 2, label = "t = $t")
end

plot!(x, y0, label = "t = 0 (exact)", color = "black", linewidth = 1, line = :dash)
for i = 1:3
    t = t0 + N*i*dt0
    y_exact = gauss.(x .- a*t)
    plot!(x, y_exact, color = "black", linewidth = 1, line = :dash, label = "t = $t (exact)")
end
display(plt)

println("\ndone")
