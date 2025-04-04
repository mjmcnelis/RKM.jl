using Revise, OrdinaryDiffEq, StaticArrays, BenchmarkTools
# using DoubleFloats: Double64
import RKM: RKM_root
using Plots;
plotly();
!(@isdefined f_ord) ? include("$RKM_root/validation/ode/logistic/equations.jl") : nothing
!(@isdefined f_ord_static) ? include("$RKM_root/validation/ode/logistic/equations.jl") : nothing
# benchmarks (alg = RK4())
# RKM  | OrdinaryDiffEq time (ms)
# 37.6 | 78.7    # save_solution = true, static_array = false
# 25.9 | 25.8    # save_solution = true, static_array = true
# 27.7 | 32.4    # save_solution = false, static_array = false
# 19.5 | 10.2    # save_solution = false, static_array = true

t0 = -10.0
N = 2
p = [0.5 - 0.25*(i-1.0)/(N-1.0+eps(1.0)) for i in 1:N]
y0 = Float64[]
for i = eachindex(p)
    push!(y0, exp(t0) / (1.0 + exp(t0)) - p[i])
end

# alg = Trapezoid(autodiff = true)
# alg = ImplicitEuler(autodiff = true)
# alg = TRBDF2(autodiff=true)
alg = RK4()

prob = ODEProblem(f_ord, y0, (t0, 10.0), p) #=f_ord_static=#
@time sol = solve(prob, alg, dt = 1e-4, reltol = 1e-6, #abstol = 0.0,
    # if adaptive is false will compute Jacobian at each timestep
    adaptive=false, saveat = 1e-4,
    # save_everystep = false, save_start = false, save_end = false
    # qmin = 0.2, qmax = 10.0, gamma = 0.9,
    # controller = IController(), dtmin = 0.0,
)

# plot(sol, legend = :outertopright) |> display
# plot!(sol.t,  mapreduce(permutedims, vcat, sol.u)) |> display

@show sol.destats
GC.gc()
println("\ndone")

# t = sol.t
# dt = t[2:end] .- t[1:end-1]
# plot(t[1:end-1], dt; size = (900, 600), linewidth = 2,
#      legend = :outertopright, legendfontsize = 12,
#      ylabel = "Δt", yguidefontsize = 14, ytickfontsize = 12,
#      xlabel = "t", xguidefontsize = 14, xtickfontsize = 12,) |> display
