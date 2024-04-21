using Revise, RKM, BenchmarkTools
import DoubleFloats: Double64
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/logistic/equations.jl") : nothing
include("$RKM_root/validation/ode/logistic/parameters.jl")

precision = Float64
# precision = Double64
# precision = BigFloat

# TODO: do asserts between adaptive, method in parameters outer-constructor
options = SolverOptions(options)

# @code_warntype TimeStepController(precision; pid = PIControl(), limiter = SmoothLimiter())
# q()

t0 = -10.0
tf = 10.0
dt0 = 1e-4

N = 2
y0 = Float64[]
for i = 1:N
    local a = N == 1 ? 0.5 : 0.5 - 0.25*(i-1.0)/(N-1.0)
    push!(y0, exp(t0) / (1.0 + exp(t0)) - a)
end

@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, precision)
# in-place version
#=
sol = Solution(precision)
@time evolve_ode!(sol, y0, t0, tf, dt0, dy_dt!, options)
=#

# @btime sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, precision)

get_stats(sol)
# plot_ode(sol, options.method, Plots.plot)

GC.gc()
println("\ndone")
