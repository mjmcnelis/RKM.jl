
using Revise
using RKM
using LinearAlgebra
using Plots
using UnPack
# plotly()

function shock(x)
    x < 0.0 ? 1.0 :
    x > 0.0 ? 0.0 : 0.5
end

# initial condition
x = LinRange(-5, 5, 101) |> collect
const dx = x[2] - x[1]
const dt = 0.05
const C = dt/(2dx)      # constant speed = (Δy/2)
@show C dx dt
N = 40

F(y) = y^2/2.0

function dy_dt!(f, t, y)
    L = length(y)
    for i in 1:L
        m = max(i-1, 1) # BC: y[0] = y[1]
        p = min(i+1, L) # BC: y[L+1] = y[L]
        ym, yc, yp = y[m], y[i], y[p]

        aR = (yc + yp)/2.0 |> abs   # |A_{i+1/2}|
        aL = (ym + yc)/2.0 |> abs   # |A_{i-1/2}|

        f[i] = -(F(yp) - F(ym))/(2.0*dx) + (aR*(yp-yc) - aL*(yc-ym))/(2.0*dx)
    end
    nothing
end

adaptive   = Fixed()
method     = Euler1()
t_range    = TimeRange(; t0 = 0.0, tf = 6.0, dt0 = dt)
options = SolverOptions(; adaptive, method, t_range)

@unpack t0, dt0 = t_range
y0 = shock.(x)

@time sol = evolve_ode(y0, dy_dt!; parameters)

plt = plot(x, y0, label = "t = 0", color = "indianred", linewidth = 2,
           size = (900, 600), ylims = (-0.25, 1.25),
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
    y_exact = shock.(x .- t/2.0)
    plot!(x, y_exact, color = "black", linewidth = 1, line = :dash, label = "t = $t (exact)")
end
# display(plt)
savefig(plt, joinpath(RKM_root, "demo/burgers_inviscid/shock/forward/upwind.png"))

println("\ndone")
