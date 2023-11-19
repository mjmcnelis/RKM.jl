
using Revise
using RKM
using LinearAlgebra
using Plots
using UnPack
# plotly()

function rarefaction(x, t)
    x < -t ? -1.0 :
    x > t ? 1.0 : x/(t + 1e-16)
end

# initial condition
x = LinRange(-5, 5, 101) |> collect
const dx = x[2] - x[1]
const dt = 0.05
const C = dt/dx         # max speed = (Δy/2)
@show C dx dt
N = 20

F(y) = y^2/2.0

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
y0 = rarefaction.(x, t0)

@time sol = evolve_ode(y0, dy_dt!; parameters)

plt = plot(x, y0, label = "t = 0", color = "indianred", linewidth = 2,
           size = (900, 600), ylims = (-1.25, 1.25),
           ylabel = "u", yguidefontsize = 14, ytickfontsize = 12,
           xlabel = "x", xguidefontsize = 14, xtickfontsize = 12,
           legend = :outertopright, legendfontsize = 12, dpi = 200)
for i = 1:3
    t = t0 + N*i*dt0
    y = sol.y[1 + N*i]
    plot!(x, y, color = "indianred", linewidth = 2, label = "t = $t")
end

plot!(x, y0, label = "t = 0 (exact)", color = "black", linewidth = 1, line = :dash)
for i = 1:3
    t = t0 + N*i*dt0
    y_exact = rarefaction.(x, t)
    plot!(x, y_exact, color = "black", linewidth = 1, line = :dash, label = "t = $t (exact)")
end
# display(plt)
savefig(plt, joinpath(RKM_root, "demo/burgers_inviscid/rarefaction/forward/lax_friedrichs.png"))

println("\ndone")
