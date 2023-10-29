using Revise, RKM, JLD2, StatsBase, Test
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/vanderpol/equations.jl") : nothing
loadpath = joinpath(RKM_root, "test/runtests/ode/answers/vanderpol_answers.jld2")

@info "Starting Vanderpol test..."

# option to reset answer keys
reset_answer_keys = false

# options to plot compare answer key / OrdinaryDiffEq
plot_compare = "ans"
# plot_compare = "diffeq"
show_plot = false

y0 = [2.0, 0.0]
t0 = 0.0
tf = 3.0e3
dt0 = 1e-4

# TODO: need to interpolate solution
parameters = Parameters(; method = TrapezoidRuleBDF2(),
                          adaptive = Doubling(; epsilon = 1e-6),
                          stage_finder = ImplicitStageFinder(; jacobian_method = ForwardJacobian()),
                        )
@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, parameters)
y, t = get_solution(sol)

# save new answer keys
if reset_answer_keys
    y_ans, t_ans = get_solution(sol)
    @save loadpath y_ans t_ans
end
@load loadpath y_ans t_ans

# plot comparison
# note: hide y2, y4 and autoscale to see y1 (RKM) vs y3 (ans/diffeq)
plt = plot_ode(sol, parameters.method, Plots.plot);
if plot_compare == "ans"
    plot!(t_ans, y_ans, color = :black, linewidth = 2, line = :dot)
elseif plot_compare == "diffeq"
    using OrdinaryDiffEq
    prob = ODEProblem(dy_dt!, y0, (t0, tf))
    @time sol = solve(prob, TRBDF2(), dt = dt0, reltol = 1e-6, abstol = 1e-6)
    plot!(sol.t, mapreduce(permutedims, vcat, sol.u),
          color = :black, linewidth = 2, line = :dash)
end
show_plot ? display(plt) : nothing

# test
for j in size(y,2)
    y_col = view(y,:,j)
    y_col_ans = view(y_ans,:,j)
    @test L1dist(y_col, y_col_ans)/mean(abs, y_col_ans) < 0.01
end

GC.gc()

@info "...done"
println("")