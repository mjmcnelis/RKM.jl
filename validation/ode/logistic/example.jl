using Revise, RKM, BenchmarkTools
import DoubleFloats: Double64
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/logistic/equations.jl") : nothing
include("$RKM_root/validation/ode/logistic/parameters.jl")

precision = Float64
# precision = Double64
# precision = BigFloat

# TODO: do asserts between adaptive, method in parameters outer-constructor
parameters = Parameters(options)

# should put t0, tf, precision here

t0 = parameters.t_range.t0

N = 2
y0 = Float64[]
for i = 1:N
    local a = N == 1 ? 0.5 : 0.5 - 0.25*(i-1.0)/(N-1.0)
    push!(y0, exp(t0) / (1.0 + exp(t0)) - a)
end

# show_progress = true
show_progress = false

static_array = false
# static_array = true

save_solution = true
# save_solution = false

@time sol = evolve_ode(y0, dy_dt!, parameters; precision, save_solution,
                                               show_progress, static_array)

# in-place version
# sol = Solution(; precision, save_solution)
# @time evolve_ode!(sol, y0, dy_dt!, parameters; show_progress, static_array)

get_stats(sol)
# plot_ode(sol, parameters.method, Plots.plot)

GC.gc()
println("\ndone")
