using Revise, RKM, OrdinaryDiffEq, LinearSolve
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/damped_harmonic_oscillator/equations.jl") : nothing

y0 = [1.0, -1.0]
t0 = 0.0 
tf = 10.0

# setup stiff ODE 
γ = 1000000000.0            # my solver doesn't work well for this large value
# γ = 1000.0                    
ω = sqrt(1.001*γ)
p = [γ, ω]

dt0 = 1e-2

# RKM
method = TrapezoidRuleBDF2()
# method = BackwardEuler1()
# method = Heun2()            # unstable
# method = Ketcheson4()       # barely stable for dt = 1e-2
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
      color = :black, linewidth = 2, line = :dash)# |> display

GC.gc()
println("\ndone")
