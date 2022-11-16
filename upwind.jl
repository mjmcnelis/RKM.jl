
using Revise
using RKM
using LinearAlgebra
using Plots 
using UnPack
plotly()

function gauss(x)
    exp(-(x - 2)^2)
end

# initial condition 
x = LinRange(0, 10, 101) |> collect
const a = 1.0
const dx = x[2] - x[1]

# upwind forward time 
function dy_dt!(f, t, y)
    f[1] = 0.0
    for i in 2:length(y)
        f[i] = -(y[i] - y[i-1]) * a/dx
    end
    nothing
end

adaptive   = Fixed()
method     = Euler1()
t_span     = TimeSpan(; t0 = 0.0, tf = 6.0, dt0 = 0.05)     # website used dt = 0.05
parameters = Parameters(; adaptive, method, t_span)

@unpack t0, dt0 = t_span
y0 = gauss.(x)
C = a*dt0/dx
@show C 

@time sol = evolve_ode(y0, dy_dt!; parameters)

plt = plot(x, y0, color = "black", size = (900, 600), linewidth = 1, line = :dash,
           ylabel = "u", yguidefontsize = 14, ytickfontsize = 12,
           xlabel = "t", xguidefontsize = 14, xtickfontsize = 12, 
           label = "t = 0 (exact)", legend = :outertopright, legendfontsize = 12)   
for n = 1:3
    t = t0 + 40*n*dt0 
    y_exact = gauss.(x .- a*t)
    plot!(x, y_exact, color = "black", linewidth = 1, line = :dash, label = "t = $t (exact)")
end

plot!(x, y0, color = "indianred", linewidth = 2, label = "t = 0 (FTBS)")
for n = 1:3
    t = t0 + 40*n*dt0 
    y = sol.y[1 + 40*n]
    plot!(x, y, color = "indianred", linewidth = 2, label = "t = $t (FTBS)")
end
display(plt)

println("\ndone")

