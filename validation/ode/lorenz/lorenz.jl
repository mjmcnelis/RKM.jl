using Revise, RKM, OrdinaryDiffEq, LinearSolve
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/lorenz/equations.jl") : nothing

y0 = [1.0, 0.0, 0.0]
t0 = 0.0
tf = 100.0

p = [10.0, 28.0, 8.0/3.0]
dt0 = 1e-3
# epsilon = 1e-8

# RKM
method = RungeKutta4()
options = SolverOptions(; method,
                          adaptive = Fixed(),
                          t_range = TimeRange(; t0, tf, dt0)
                        )
@time sol = evolve_ode(y0, dy_dt!; model_parameters = p, parameters, show_progress = false)
get_stats(sol)
plot_ode(sol, method, Plots.plot);

# OrdinaryDiffEq
prob = ODEProblem(dy_dt!, y0, (t0, tf), p)
@time sol = solve(prob, RK4(), reltol = epsilon, abstol = epsilon, dt = dt0, adaptive = false)
@show sol.destats
plot!(sol.t,  mapreduce(permutedims, vcat, sol.u),
      color = :black, linewidth = 2, line = :dash) |> display

GC.gc()
println("\ndone")
