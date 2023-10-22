using Revise
using RKM, Test
# using Plots; plotly()
import DoubleFloats: Double64
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/logistic/equations.jl") : nothing

precision_vect = [Float32, Float64, Double64]   #
method_vect = [HeunEuler21()]
adaptive_vect = [Fixed(), Doubling(), Embedded()]

t0 = -8.0
tf = 8.0
dt0 = 1e-4

y0 = exp(t0)/(1.0 + exp(t0)) - 0.5

for precision in precision_vect
    for method in method_vect
        for adaptive in adaptive_vect
            save_solution = adaptive isa Fixed
            # save_solution = true

            # @show precision adaptive

            parameters = Parameters(; method, adaptive, save_solution)

            # sol = Solution(; precision)
            # evolve_ode!(sol, y0, t0, tf, dt0, dy_dt!, parameters)
            # @show sol.excess_allocations[1]
            # note: second solve does not allocate on adaptive even if save_solution = true
            # evolve_ode!(sol, y0, t0, tf, dt0, dy_dt!, parameters)

            sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, parameters)
            # @show sol.excess_allocations[1]
            sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, parameters)

            @test sol.excess_allocations[1] == 0
            # @show sol.excess_allocations[1]

            # note: plotting here gives me this error unless deepcopy
            # ERROR: LoadError: cannot resize array with shared data
            # plot_ode(deepcopy(sol), parameters.method, Plots.plot) |> display
            # println("")
        end
    end
end

println("\ndone")
