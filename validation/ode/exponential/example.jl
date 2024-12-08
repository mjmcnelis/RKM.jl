using Revise, RKM, BenchmarkTools
import DoubleFloats: Double64
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/exponential/equations.jl") : nothing

# p = [1.0]           # exponential growth
p = [-1.0]          # exponential decay

t0 = 0.0 |> BigFloat
y0 = exp(p[1] * t0)

tf = 10.0
dt0 = 1e-4

options = SolverOptions(method = RungeKutta4(),
                        adaptive = Fixed(),
                        sensitivity_method = DecoupledDirect()
                       )

@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p)

get_stats(sol)

plot_ode(sol, options.method, Plots.plot)#; logy = true, show_time_step = true)
plot!(sol.t, y_exact.(sol.t; p), label = "exact") |> display

t, S = get_sensitivity(sol)
plot(t, S) |> display

GC.gc()
println("\ndone")
