using Revise, RKM, BenchmarkTools, UnPack, Setfield
import DoubleFloats: Double64
using LinearSolve, RecursiveFactorization, ProgressMeter
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

# for sensitivity testing only
# sparse_Jy = nansafe_state_jacobian(y0, t0, dy_dt!, p; chunk_size = 1);
# sparse_Jp = nansafe_param_jacobian(y0, t0, dy_dt!, p; chunk_size = 1);
# @set! options.stage_finder.state_jacobian = ForwardColorJacobian(; sparsity = sparse_Jy)
# @set! options.sensitivity.param_jacobian = ForwardColorJacobian(; sparsity = sparse_Jp)

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
#=
plot(t_dense, y_dense; color = [:red :blue], linewidth = 2, legend = :outertopright);
plot!(t_dense, hcat(y_exact.(t_dense; N = 2)...)';
      color = :black, linewidth = 2, line = :dash)
=#

GC.gc()
println("\ndone")
