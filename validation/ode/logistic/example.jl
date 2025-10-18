using Revise, RKM, BenchmarkTools, UnPack, Setfield
import DoubleFloats: Double64
using LinearSolve, RecursiveFactorization
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/logistic/equations.jl") : nothing
include("$RKM_root/validation/ode/logistic/parameters.jl")

show_plot = false

# TODO: do asserts between adaptive, method in parameters outer-constructor
options = SolverOptions(options)

t0 = -10.0
tf = 10.0
dt0 = 1e-4

N = 2
p = [0.5 - 0.25*(i-1.0)/(N-1.0+eps(1.0)) for i in 1:N]
y0 = Float64[]
for i = eachindex(p)
    push!(y0, exp(t0) / (1.0 + exp(t0)) - p[i])
end

y0 = Complex.(y0)

# for sensitivity testing only
# sparse_Jy = nansafe_state_jacobian(y0, t0, dy_dt!, p; chunk_size = 1);
# sparse_Jp = nansafe_param_jacobian(y0, t0, dy_dt!, p; chunk_size = 1);
# @set! options.state_jacobian = ForwardJacobian(; sparsity = sparse_Jy)
# @set! options.sensitivity.param_jacobian = ForwardJacobian(; sparsity = sparse_Jp)

@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p)
# in-place version
#=
sol = Solution(options)
@time evolve_ode!(sol, y0, t0, tf, dt0, dy_dt!, options, p)
=#

# @btime sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p)

get_stats(sol)
# plot_ode(sol, options.method, Plots.plot)

# analyze jacobian spectrum
# TODO: interpolation should follow same arguments as get_eigenvalues
@unpack method, eigenmax = options
if method.iteration isa RKM.Implicit && !(eigenmax isa NoEigenMax)
    @time t, lambda = get_eigenvalues(sol, dy_dt!, options, p; dt_dense = 1e-2)
    _, lambda_LR = get_eigenmax(sol)
    if show_plot
        plot(t, real.(lambda))
        plot!(sol.t, real.(lambda_LR), line = :dash) |> display
    end
end

# interpolation
@unpack interpolator = options
if !(interpolator isa NoInterpolation)
    @time t_dense, y_dense = interpolate_solution(options, sol; dt_dense = 1e-4)
    if show_plot
        t, y = get_solution(sol)
        scatter(t, y; color = [:red :blue], ms = 3,
                size = (900,600), legend = :outertopright)
        plot!(t_dense, y_dense; color = [:red :blue], linewidth = 2) |> display
    end
end

GC.gc()
println("\ndone")
