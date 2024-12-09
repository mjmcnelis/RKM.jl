using Revise, RKM, BenchmarkTools
import DoubleFloats: Double64
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/power/equations.jl") : nothing

# p = [2.0]           # power growth
p = [-2.0]          # power decay

t0 = 1.0 |> BigFloat
y0 = t0^p[1]

tf = 10.0
dt0 = 1e-4

options = SolverOptions(method = RungeKutta4(),
                        adaptive = Fixed(),
                        sensitivity_method = DecoupledDirect()
                       )

@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p)

get_stats(sol)

plot_ode(sol, options.method, Plots.plot)
plot!(sol.t, y_exact.(sol.t; p), label = "exact") |> display

t, S = get_sensitivity(sol)
plot(t, S) |> display

GC.gc()
println("\ndone")
