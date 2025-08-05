using Revise, RKM, OrdinaryDiffEq, LinearSolve
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/lotka_vorterra/equations.jl") : nothing

y0 = [1.0, 1.0]
t0 = 0.0
tf = 10.0

p = [1.5, 1.0, 3.0, 1.0]
dt0 = 1e-4
epsilon = 1e-8

# RKM
method = DormandPrince5()
options = SolverOptions(; method,
                        #   adaptive = Embedded(; epsilon),
                          adaptive = Fixed(),
                          sensitivity = DecoupledDirect(),
                        )
@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p)
get_stats(sol)
plot(get_sensitivity(sol)) |> display

plt = plot_ode(sol, method, Plots.plot);

# OrdinaryDiffEq
prob = ODEProblem(dy_dt!, y0, (t0, tf), p)
@time sol = solve(prob, DP5(), reltol = epsilon, abstol = epsilon, dt = dt0)
sol.destats |> display
plot!(sol.t,  mapreduce(permutedims, vcat, sol.u),
      color = :black, linewidth = 2, line = :dash)
display(plt)
GC.gc()
println("\ndone")
