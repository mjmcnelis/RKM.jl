using RKM, JLD2, StatsBase, Test
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/flame/equations.jl") : nothing
loadpath = joinpath(RKM_root, "test/ode/answers/flame_answers.jld2")
@info "Starting flame test..."

# option to reset answer keys
reset_answer_keys = false

# options to plot, compare answer key / OrdinaryDiffEq
show_plot = false
plot_compare = "ans"
# plot_compare = "diffeq"

# plot eigenvalues of state jacobian
plot_eigvals = false

δ = 1e-4
y0 = [δ]
t0 = 0.0
tf = 2.0/δ
dt0 = δ

# TODO: need to interpolate solution
options = SolverOptions(; method = TrapezoidRuleBDF2(),
                          adaptive = Doubling(; epsilon = 1e-7, alpha = 1e-7, delta = 1e-7),
                          stage_finder = ImplicitStageFinder(;
                                             state_jacobian = ForwardJacobian(),
                                         ),
                          eigenmax = KrylovEigenMax(; krylovdim = 1)
                        )
@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options)
t, y = get_solution(sol)

# might want to create a plot function to compare eigenvalues
if plot_eigvals
    using Plots; plotly()
    t, lambda = get_eigenvalues(sol, dy_dt!, options; dt_dense = 1e-2)
    plot(t, real.(lambda))
    _, lambda_LR = get_eigenmax(sol)
    plot!(sol.t, real.(lambda_LR), line = :dash) |> display
end

# save new answer keys
if reset_answer_keys
    t_ans, y_ans = get_solution(sol)
    @save loadpath t_ans y_ans
end
@load loadpath t_ans y_ans

# plot comparison
if show_plot
    using Plots; plotly()
    plt = plot_ode(sol, options.method, Plots.plot);
    if plot_compare == "ans"
        plot!(t_ans, y_ans, color = :black, linewidth = 2, line = :dot)
    elseif plot_compare == "diffeq"
        using OrdinaryDiffEq
        prob = ODEProblem(dy_dt!, y0, (t0, tf))
        @time sol = solve(prob, TRBDF2(), dt = dt0, reltol = 1e-7, abstol = 1e-7)
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