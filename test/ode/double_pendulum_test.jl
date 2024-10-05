using RKM, JLD2, StatsBase, Test
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/double_pendulum/equations.jl") : nothing
loadpath = joinpath(RKM_root, "test/ode/answers/double_pendulum_answers.jld2")
@info "Starting double pendulum test..."

# option to reset answer keys
reset_answer_keys = false

# options to plot, compare answer key / OrdinaryDiffEq
show_plot = false
plot_compare = "ans"
# plot_compare = "diffeq"

y0 = [pi/2.0, pi/2.0, 0.0, 0.0]
t0 = 0.0
tf = 10.0

p = [0.5]
dt0 = 1e-4

options = SolverOptions(; method = RungeKutta4(), adaptive = Fixed())
@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p)
y, t = get_solution(sol)

# save new answer keys
if reset_answer_keys
    y_ans, t_ans = get_solution(sol)
    @save loadpath y_ans t_ans
end
@load loadpath y_ans t_ans

# plot comparison
if show_plot
    using Plots; plotly()
    plt = plot_ode(sol, options.method, Plots.plot);
    if plot_compare == "ans"
        plot!(t_ans, y_ans, color = :black, linewidth = 2, line = :dot)
    elseif plot_compare == "diffeq"
        using OrdinaryDiffEq
        prob = ODEProblem(dy_dt!, y0, (t0, tf), p)
        @time sol = solve(prob, RK4(), dt = dt0, adaptive = false)
        plot!(sol.t, mapreduce(permutedims, vcat, sol.u),
              color = :black, linewidth = 2, line = :dash)
    end
    display(plt)
end

# test
for j in size(y,2)
    y_col = view(y,:,j)
    y_col_ans = view(y_ans,:,j)
    @test L1dist(y_col, y_col_ans)/mean(abs, y_col_ans) < 0.01
end

GC.gc()

@info "...done"
println("")