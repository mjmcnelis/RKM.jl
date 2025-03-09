using Revise, RKM, BenchmarkTools, UnPack
import DoubleFloats: Double64
using LinearSolve, RecursiveFactorization
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/logistic/equations.jl") : nothing
include("$RKM_root/validation/ode/logistic/parameters.jl")

# TODO: do asserts between adaptive, method in parameters outer-constructor
options = SolverOptions(options)

# @code_warntype TimeStepController(precision; pid = PIControl(), limiter = SmoothLimiter())
# q()

t0 = -10.0
tf = 10.0
dt0 = 1e-4

N = 2
p = [0.5 - 0.25*(i-1.0)/(N-1.0+eps(1.0)) for i in 1:N]
y0 = Float64[]
for i = eachindex(p)
    push!(y0, exp(t0) / (1.0 + exp(t0)) - p[i])
end

@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p)
# in-place version
#=
sol = Solution(options)
@time evolve_ode!(sol, y0, t0, tf, dt0, dy_dt!, options)
=#

# @btime sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options)

get_stats(sol)
# plot_ode(sol, options.method, Plots.plot)

# interpolation
@time t_dense, y_dense = interpolate_solution(options, sol; dt_dense = 1e-5)

GC.gc()
println("\ndone")
