
using Revise
using RKM
using LinearAlgebra
using Plots
using UnPack
plotly()

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
        f[i] = -(F(yp) - F(ym))/(2.0*dx) + (yp - 2.0*yc + ym)/(2.0*dt)
    end
    nothing
end

function jacobian!(J, t, y)
    nrow = size(J,1)
    A = 1.0/(2.0*dx)
    B = 1.0/(2.0*dt)
    # note: includes NBC
    J[1,1]       = A*y[1] - B
    J[1,2]       = B - A*y[2]
    J[end,end]   = -A*y[end] - B
    J[end,end-1] = A*y[end-1] + B
    for i in 2:nrow-1
        J[i,i-1] = A*y[i-1] + B
        J[i,i]   = -2.0*B
        J[i,i+1] = B - A*y[i+1]
    end
    nothing
end

adaptive   = Fixed()
method     = BackwardEuler1()
t_range    = TimeRange(; t0 = 0.0, tf = 6.0, dt0 = dt)
parameters = Parameters(; adaptive, method, t_range)

@unpack t0, dt0 = t_range
y0 = shock.(x)

@time sol = evolve_ode(y0, dy_dt!; jacobian!, parameters)

plt = plot(x, y0, label = "t = 0", color = "indianred", linewidth = 2,
           size = (900, 600), ylims = (-0.25, 1.25),
           ylabel = "u", yguidefontsize = 14, ytickfontsize = 12,
           xlabel = "x", xguidefontsize = 14, xtickfontsize = 12,
           legend = :outertopright, legendfontsize = 12)
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
display(plt)

println("\ndone")
