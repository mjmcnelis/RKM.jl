using Revise, RKM, OrdinaryDiffEq
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/vanderpol/equations.jl") : nothing

y0 = [2.0, 0.0]
t0 = 0.0 
tf = 3.0e3 
dt0 = 1e-4
epsilon = 1e-6

parameters = Parameters(; method = TrapezoidRuleBDF2(), 
                          adaptive = Doubling(; epsilon), 
                          stage_finder = ImplicitStageFinder(; jacobian_method = ForwardJacobian()),
                          t_range = TimeRange(; t0, tf, dt0)
                        )
@time sol = evolve_ode(y0, dy_dt!; parameters)
get_stats(sol)

# note: hide y2, y4 and autoscale to see y1 vs y3
plot_ode(sol, method, Plots.plot);

epsilon = 1e-7
prob = ODEProblem(dy_dt!, y0, (t0, tf))
@time sol = solve(prob, TRBDF2(), dt = dt0, reltol = epsilon, abstol = epsilon)

if sol.retcode == ReturnCode.Success
    display(plot!(sol.t,  mapreduce(permutedims, vcat, sol.u), 
                  color = :black, linewidth = 2, line = :dash)
            )
    @show sol.destats
else 
    @show sol.retcode
end

GC.gc()
println("\ndone")
