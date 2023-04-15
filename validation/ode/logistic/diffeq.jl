using Revise, OrdinaryDiffEq, StaticArrays, BenchmarkTools
using DoubleFloats: Double64
import RKM: RKM_root
using Plots; plotly()
!(@isdefined f_ord) ? include("$RKM_root/validation/ode/logistic/equations.jl") : nothing 

t0 = -10.0
N = 1 
y0 = [exp(t0) / (1.0 + exp(t0)) - get_a(i,N) for i = 1:N]

prob = ODEProblem(f_ord, y0, (t0, 10.0))
@time sol = solve(prob, Heun(), dt = 1e-4, reltol = 1e-6, abstol = 0.0, 
                  # qmin = 0.2, qmax = 10.0, gamma = 0.9,
                  # controller = IController(), dtmin = 0.0,
                  )

plot(sol, legend = :outertopright) |> display

@show sol.destats
GC.gc()
println("\ndone")

t = sol.t
dt = t[2:end] .- t[1:end-1]
plot(t[1:end-1], dt; size = (900, 600), linewidth = 2,
     legend = :outertopright, legendfontsize = 12,
     ylabel = "Δt", yguidefontsize = 14, ytickfontsize = 12,
     xlabel = "t", xguidefontsize = 14, xtickfontsize = 12,) |> display
     