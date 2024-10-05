using RKM, Test
# using Plots; plotly()
import DoubleFloats: Double64
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/logistic/equations.jl") : nothing

@info "Starting allocation test..."

precision_vect = [Float64, Double64]
method_vect = [HeunEuler21()]
adaptive_vect = [Fixed(), Doubling(), Embedded()]
limiter_vect = [PiecewiseLimiter(), SmoothLimiter()]
pid_vect = [BasicControl(), PIControl(), H312Control(),
            H321PredictiveControl(), H211bPredictiveControl()]

t0 = -8.0
tf = 8.0
dt0 = 1e-4

p = [0.5]
y0 = exp(t0)/(1.0 + exp(t0)) - p[1]

for precision in precision_vect, method in method_vect,
                                 adaptive in adaptive_vect,
                                 limiter in limiter_vect,
                                 pid in pid_vect

    controller = TimeStepController(precision; pid, limiter)
    save_solution = adaptive isa Fixed
    # save_solution = true

    # @show precision adaptive limiter pid

    options = SolverOptions(; method, adaptive, controller, save_solution, precision)

    # sol = Solution(precision)
    # evolve_ode!(sol, y0, t0, tf, dt0, dy_dt!, options)
    # @show sol.excess_memory[1]
    # note: second solve does not allocate on adaptive even if save_solution = true
    # evolve_ode!(sol, y0, t0, tf, dt0, dy_dt!, options)

    sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p)
    # @show sol.excess_memory[1]
    sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p)

    @test sol.excess_memory[1] == 0
    # @show sol.excess_memory[1]

    # note: plotting here gives me this error unless deepcopy
    # ERROR: LoadError: cannot resize array with shared data
    # plot_ode(deepcopy(sol), options.method, Plots.plot) |> display
    # println("")
end

@info "...done"
println("")
