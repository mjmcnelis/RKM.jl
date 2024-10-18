using RKM, JLD2, StatsBase, Test
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/lorenz/equations.jl") : nothing
loadpath = joinpath(RKM_root, "test/ode/answers/lorenz_answers.jld2")
@info "Starting Lorenz test..."

# option to reset answer keys
reset_answer_keys = false

# options to plot, compare answer key / OrdinaryDiffEq
show_plot = false
plot_compare = "ans"
# plot_compare = "diffeq"

if show_plot && plot_compare == "diffeq"
    precision = BigFloat
else
    precision = Float64
end

y0 = [1.0, 0.0, 0.0]
t0 = 0.0
tf = 100.0

p = [10.0, 28.0, precision(8//3)]
dt0 = 1e-3

options = SolverOptions(; method = Feagin108(), adaptive = Fixed(), precision)
@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p)
t, y = get_solution(sol)

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
        y0 = y0 .|> precision
        t0 = t0 |> precision
        tf = tf |> precision
        dt0 = dt0 |> precision
        prob = ODEProblem(dy_dt!, y0, (t0, tf), p)
        @time sol = solve(prob, Feagin10(), dt = dt0, adaptive = false)
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