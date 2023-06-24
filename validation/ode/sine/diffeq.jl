using Revise, OrdinaryDiffEq, StaticArrays, BenchmarkTools
import RKM: RKM_root
# using Plots; plotly()
!(@isdefined fp) ? include("$RKM_root/validation/ode/sine/equations.jl") : nothing

t0 = 0.0
u0 = [sin(t0), ω*cos(ω*t0)]

prob = ODEProblem(fp, u0, (t0, pi/2))
@time sol = solve(prob, RK4(), dt = 1e-5, adaptive = false)
# @btime sol = solve(prob, RK4(), dt = 1e-5, adaptive = false)
# plot(sol, legend = :outertopright)

GC.gc()
println("\ndone")
