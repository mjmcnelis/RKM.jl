using Revise, RKM, BenchmarkTools
import DoubleFloats: Double64
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/projectile/equations.jl") : nothing

γ = 0.1                 # damping coefficient
g = 9.81                # gravitational acceleration

p = [γ, g]

t0 = 0.0
y0 = [0.0, 100.0]       # y, vy

tf = 10.0
dt0 = 1e-4

options = SolverOptions(method = RungeKutta4(), adaptive = Fixed())

@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p)
t, _ = get_solution(sol)
get_stats(sol)

plot_ode(sol, options.method, Plots.plot)
plot!(t, hcat(y_exact.(t; y0, t0, p)...)', label = "exact", linewidth = 1.5) |> display

GC.gc()
println("\ndone")
