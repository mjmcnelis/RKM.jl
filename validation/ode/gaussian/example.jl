using Revise, RKM, BenchmarkTools
import DoubleFloats: Double64
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/gaussian/equations.jl") : nothing

# precision = Double64
precision = Float64

t0 = -5.0
y0 = exp(-t0^2/2.0) + 1.0     # note: shift by +1 to ignore abstol

adaptive = Embedded(; epsilon = 1e-15, low = 0.1, high = 5.0, safety = 0.9)
method = DormandPrince54()
t_range = TimeRange(; t0, tf = 5, dt0 = 1e-4)
timer = TimeLimit()
options = SolverOptions(; adaptive, method, t_range, timer)

@time sol = evolve_ode(y0, dy_dt!; options, precision, show_progress = false)
# plot_ode(sol, method, Plots.plot)

@show timer.counter sol.FE sol.rejection_rate
GC.gc()
println("\ndone")

y, t = get_solution(sol)

dt = t[2:end] .- t[1:end-1]
plot!(t[1:end-1], dt; size = (900, 600), linewidth = 2,
     legend = :outertopright, legendfontsize = 12,
     ylabel = "Δt", yguidefontsize = 14, ytickfontsize = 12,
     xlabel = "t", xguidefontsize = 14, xtickfontsize = 12,) |> display
