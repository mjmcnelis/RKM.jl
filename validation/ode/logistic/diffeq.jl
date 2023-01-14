using Revise, OrdinaryDiffEq, StaticArrays, BenchmarkTools
import RKM: RKM_root
# using Plots; plotly()
!(@isdefined fp) ? include("$RKM_root/validation/ode/logistic/equations.jl") : nothing

t0 = -10.0
N = 2
# u0 = SA[[exp(t0)/(1.0 + exp(t0)) - get_a(i,N) for i = 1:N]...]
u0 = Float64[[exp(t0)/(1.0 + exp(t0)) - get_a(i,N) for i = 1:N]...]

prob = ODEProblem(fp, u0, (t0, 10.0))
@time sol = solve(prob, RK4(), dt = 1e-4, adaptive = false)
# @btime sol = solve(prob, RK4(), dt = 1e-4, adaptive = false)
# plot(sol, legend = :outertopright)

GC.gc()
println("\ndone")
