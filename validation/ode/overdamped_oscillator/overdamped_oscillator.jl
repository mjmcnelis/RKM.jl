using Revise, RKM, OrdinaryDiffEq, LinearSolve
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/overdamped_oscillator/equations.jl") : nothing

y0 = [1.0, -1.0]        # eigenvector of ODE system (exact solution is y(t) = exp(-t)*y0)
t0 = 0.0
tf = 10.0

# setup stiff ODE
# ω = 10.0^(9/2)          # my solver is slow w/ this large value when epsilon = 1e-8
ω = 10.0^(3/2)
γ = ω^2 + 1.0           # eigenvalues of ODE are λ = (-ω², -1)

p = [γ, ω]

dt0 = 1e-2

# RKM
method = TrapezoidRuleBDF2()
# method = BackwardEuler1()
# method = Heun2()            # unstable
# method = Ketcheson4()       # barely stable for ω² = 1000, dt = 0.01
jacobian_method = ForwardJacobian()

parameters = Parameters(; method,
                          adaptive = Fixed(),
                          stage_finder = ImplicitStageFinder(; epsilon = 1e-6, # get better performance
                                                               jacobian_method,
                                                            ),
                          t_range = TimeRange(; t0, tf, dt0)
                        )
@time sol = evolve_ode(y0, dy_dt!; model_parameters = p, parameters, show_progress = false)
get_stats(sol)
plot_ode(sol, method, Plots.plot);

# OrdinaryDiffEq
prob = ODEProblem(dy_dt!, y0, (t0, tf), p)
@time sol = solve(prob, TRBDF2(linsolve = LUFactorization()), dt = dt0, adaptive = false)
@show sol.destats
plot!(sol.t,  mapreduce(permutedims, vcat, sol.u),
      color = :black, linewidth = 2, line = :dash) |> display

GC.gc()
println("\ndone")
