# note: used to benchmark performance against OrdinaryDiffEq
using Revise
using OrdinaryDiffEq
import RKM: RKM_root
using StaticArrays
using BenchmarkTools
using Plots
plotly()
!(@isdefined fp) ? include("$RKM_root/equations.jl") : nothing 

t0 = -10.0
N = 2

u0 = SA[[exp(t0)/(1.0 + exp(t0)) - get_a(i,N) for i = 1:N]...]
# u0 = Float64[[exp(t0)/(1.0 + exp(t0)) - get_a(i,N) for i = 1:N]...]

# prob = ODEProblem(fp, u0, (t0, t0))
# integ = init(prob, RK4(), dt = 1e-2, adaptive = false)#, saveat = 1000)
# @time solve!(integ)

GC.gc()
prob = ODEProblem(fp, u0, (t0, 10.0))
@btime sol = solve(prob, RK4(), dt = 1e-4, adaptive = false)#, saveat = 1000)
# @show sol.destats
println("\ndone")

# @benchmark sol = solve(prob, RK4(), dt = 1e-4, adaptive = false)

#=
plot(sol,linewidth=5,title="Solution to the linear ODE with a thick line",
     xaxis="Time (t)",yaxis="u(t) (in μm)",label="My Thick Line!") # legend=false
=#
