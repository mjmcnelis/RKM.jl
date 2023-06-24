using Revise, RKM, OrdinaryDiffEq
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/robertson/equations.jl") : nothing

y0 = [1.0, 0.0, 0.0]
t0 = 0.01
tf = 1.0e4
dt0 = 0.01

epsilon = 1e-6

# RKM
method = TrapezoidRuleBDF2()
parameters = Parameters(; method,
                          adaptive = Doubling(; epsilon),
                          stage_finder = ImplicitStageFinder(; jacobian_method = ForwardJacobian()),
                          t_range = TimeRange(; t0, tf, dt0)
                        )
@time sol = evolve_ode(y0, dy_dt!; parameters)
get_stats(sol)
plot_ode(sol, method, Plots.plot; logx = true)

# OrdinaryDiffEq
prob = ODEProblem(dy_dt!, y0, (t0, tf))
@time sol = solve(prob, TRBDF2(), dt = dt0, reltol = epsilon, abstol = epsilon)
@show sol.destats
plot!(sol.t,  mapreduce(permutedims, vcat, sol.u),
      color = :black, linewidth = 2, line = :dash) |> display

# diffeq goes her


GC.gc()
println("\ndone")
