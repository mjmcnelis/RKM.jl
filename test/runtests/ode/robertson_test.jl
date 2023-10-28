using Revise, RKM, JLD2, StatsBase, Test
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/robertson/equations.jl") : nothing
loadpath = joinpath(RKM_root, "test/runtests/ode/roberston_answers.jld2")

@info "Starting Robertson test..."

# option to reset answer keys
reset_answer_keys = false

# options to plot compare answer key / OrdinaryDiffEq
plot_compare = "ans"
# plot_compare = "diffeq"
show_plot = false

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

@load loadpath y_ans t_ans

# plot comparison
plt = plot_ode(sol, parameters.method, Plots.plot; logx = true);
if plot_compare == "ans"
    plot!(t_ans, y_ans, color = :black, linewidth = 2, line = :dot)
elseif plot_compare == "diffeq"
    using OrdinaryDiffEq
    prob = ODEProblem(dy_dt!, y0, (t0, tf))
    @time sol = solve(prob, TRBDF2(), dt = dt0, reltol = 1e-6, abstol = 1e-6);
    plot!(sol.t,  mapreduce(permutedims, vcat, sol.u),
          color = :black, linewidth = 2, line = :dash)
end
show_plot ? display(plt) : nothing

# test
for j in size(y,2)
    y_col = view(y,:,j)
    y_col_ans = view(y_ans,:,j)
    @test L1dist(y_col, y_col_ans)/mean(abs, y_col_ans) < 0.01
end

# save new answer keys
if reset_answer_keys
    y_ans, t_ans = get_solution(sol)
    @save loadpath y_ans t_ans
end

GC.gc()

@info "...done"
println("")