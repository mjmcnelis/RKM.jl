using Revise, RKM, JLD2, StatsBase, Test
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/robertson/equations.jl") : nothing
loadpath = joinpath(RKM_root, "test/runtests/ode/roberston_answers.jld2")

y0 = [1.0, 0.0, 0.0]
t0 = 0.01
tf = 1.0e4
dt0 = 0.01

# TODO: need to interpolate solution
parameters = Parameters(; method = TrapezoidRuleBDF2(),
                          adaptive = Doubling(; epsilon = 1e-6),
                          stage_finder = ImplicitStageFinder(; jacobian_method = ForwardJacobian()),
                        )
@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, parameters)
y, t = get_solution(sol)

#=
y_ans, t_ans = get_solution(sol)
@save loadpath y_ans t_ans
=#
@load loadpath y_ans t_ans

for j in size(y,2)
    y_col = view(y,:,j)
    y_col_ans = view(y_ans,:,j)
    @test L1dist(y_col, y_col_ans)/mean(abs, y_col_ans) < 0.01
end

plot_ode(sol, parameters.method, Plots.plot; logx = true)

# for comparison to OrdinaryDiffEq
#=

using OrdinaryDiffEq
prob = ODEProblem(dy_dt!, y0, (t0, tf))
@time sol = solve(prob, TRBDF2(), dt = dt0, reltol = 1e-6, abstol = 1e-6)
@show sol.destats
plot!(sol.t,  mapreduce(permutedims, vcat, sol.u),
      color = :black, linewidth = 2, line = :dash) |> display
=#
GC.gc()
println("\ndone")
