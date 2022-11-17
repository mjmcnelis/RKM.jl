
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
N = 40

# lax-friedrichs forward time 
function dy_dt!(f, t, y)
    dt = 0.05 # TEMP 
    
    # y[0] = y[1]
    fR = (y[1] + y[2])*a/2.0 - (y[2] - y[1])*dx/(2.0*dt)
    fL = y[1]*a
    f[1] = -(fR - fL) / dx
    for i in 2:length(y)-1
        fR = (y[i] + y[i+1])*a/2.0 - (y[i+1] - y[i])*dx/(2.0*dt)
        fL = (y[i-1] + y[i])*a/2.0 - (y[i] - y[i-1])*dx/(2.0*dt)
        f[i] = -(fR - fL) / dx
    end
    # y[end+1] = y[end]
    fR = y[end]*a
    fL = (y[end-1] + y[end])*a/2.0 - (y[end] - y[end-1])*dx/(2.0*dt)
    f[end] = -(fR - fL) / dx
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

plt = plot(x, y0, label = "t = 0", color = "indianred", linewidth = 2,
           size = (900, 600), ylims = (-0.5, 1.3),
           ylabel = "u", yguidefontsize = 14, ytickfontsize = 12, 
           xlabel = "t", xguidefontsize = 14, xtickfontsize = 12, 
           legend = :outertopright, legendfontsize = 12)
for i = 1:3
    t = t0 + N*i*dt0 
    y = sol.y[1 + N*i]
    plot!(x, y, color = "indianred", linewidth = 2, label = "t = $t")
end

plot!(x, y0, label = "t = 0 (exact)", color = "black", linewidth = 1, line = :dash)   
for i = 1:3
    t = t0 + N*i*dt0 
    y_exact = gauss.(x .- a*t)
    plot!(x, y_exact, color = "black", linewidth = 1, line = :dash, label = "t = $t (exact)")
end
display(plt)

println("\ndone")

