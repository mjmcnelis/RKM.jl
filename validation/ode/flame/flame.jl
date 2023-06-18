using Revise, RKM, OrdinaryDiffEq
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/flame/equations.jl") : nothing

δ = 1e-4
y0 = [δ]
t0 = 0.0 
tf = 2.0/δ
dt0 = δ

epsilon = 1e-7

# RKM
method = TrapezoidRuleBDF2()
parameters = Parameters(; method, 
                          adaptive = Doubling(; epsilon), 
                          stage_finder = ImplicitStageFinder(; jacobian_method = ForwardJacobian()),
                          t_range = TimeRange(; t0, tf, dt0)
                        )
@time sol = evolve_ode(y0, dy_dt!; parameters)
get_stats(sol)
plot_ode(sol, method, Plots.plot);

# OrdinaryDiffEq
prob = ODEProblem(dy_dt!, y0, (t0, tf))
@time sol = solve(prob, TRBDF2(), dt = dt0, reltol = epsilon, abstol = epsilon)
@show sol.destats
plot!(sol.t,  mapreduce(permutedims, vcat, sol.u), 
      color = :black, linewidth = 2, line = :dash) |> display

GC.gc()
println("\ndone")
