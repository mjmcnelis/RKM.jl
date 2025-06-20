using RKM, JLD2, StatsBase, Test
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/overdamped_oscillator/equations.jl") : nothing
loadpath = joinpath(RKM_root, "test/ode/answers/overdamped_oscillator_answers.jld2")
@info "Starting overdamped oscillator test..."

# option to reset answer keys
reset_answer_keys = false

# options to plot, compare answer key / OrdinaryDiffEq
show_plot = false
plot_compare = "ans"
# plot_compare = "diffeq"

plot_eigvals = false

y0 = [1.0, -1.0]    # eigenvector of ODE system (exact solution is y(t) = ±exp(-t)*y0)
t0 = 0.0
tf = 10.0
dt0 = 1e-2

# note: RKM slow for very large ω when epsilon = 1e-8 (better performance if use 1e-6)
# note: Double64 allows for 25/2, 71/2 for BigFloat
ω = 10.0^(9/2)
# ω = 10.0^(3/2)
γ = ω^2 + 1.0       # ODE eigenvalues are λ = (-ω², -1), making it very stiff
p = [γ, ω]

# TODO: need to interpolate solution
options = SolverOptions(; adaptive = Fixed(),
                          method = TrapezoidRuleBDF2(),
                        #   method = Ketcheson4(), # barely stable for ω = 10.0^(3/2)
                          state_jacobian = ForwardJacobian(),
                          root_finder = Newton(; epsilon = 1e-6),
                          eigenmax = KrylovEigenMax(; krylovdim = 2),
                       )
@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p)
t, y = get_solution(sol)

if plot_eigvals
    using Plots; plotly()
    t, lambda = get_eigenvalues(sol, dy_dt!, options, p; dt_dense = 1e-2)
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
        using OrdinaryDiffEq, LinearSolve
        prob = ODEProblem(dy_dt!, y0, (t0, tf), p)
        @time sol = solve(prob, TRBDF2(linsolve = LUFactorization()),
                        dt = dt0, adaptive = false)
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