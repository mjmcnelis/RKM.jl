using Revise, RKM, OrdinaryDiffEq, LinearSolve
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/double_pendulum/equations.jl") : nothing

y0 = [pi/2.0, pi/2.0, 0.0, 0.0]
t0 = 0.0
tf = 10.0

p = [0.5]
dt0 = 1e-4

# RKM
method = RungeKutta4()
options = SolverOptions(; method,
                          adaptive = Fixed(),
                          t_range = TimeRange(; t0, tf, dt0)
                        )
@time sol = evolve_ode(y0, dy_dt!; model_parameters = p, parameters,
                                   show_progress = false
                      )
get_stats(sol)
plot_ode(sol, method, Plots.plot);

# OrdinaryDiffEq
prob = ODEProblem(dy_dt!, y0, (t0, tf), p)
@time sol = solve(prob, RK4(), dt = dt0, adaptive = false)
@show sol.destats
plot!(sol.t,  mapreduce(permutedims, vcat, sol.u),
      color = :black, linewidth = 2, line = :dash) |> display

GC.gc()
println("\ndone")
