using Revise
using RKM
using BenchmarkTools
import DoubleFloats: Double64
import StaticArrays: SA
using Plots
plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/equations.jl") : nothing 

precision = Float64
# precision = Double64

adaptive = Fixed()             # 25 allocs (dt = 1e-4)
# adaptive = Doubling()          # 30 allocs (Heun2)
# adaptive = Embedded()          # 27 allocs (Fehlberg45)

method = RungeKutta4()
# method = Heun2()
# method = Fehlberg45()
# method = BackwardEuler1()

t_range = TimeRange(; t0 = -10.0, tf = 10.0, dt0 = 1e-4)
timer = TimeLimit()

# do asserts between adaptive, method in parameters outer-constructor
parameters = Parameters(; adaptive, method, t_range, timer)

t0 = t_range.t0

N = 2
y0 = Float64[]
for i = 1:N
    local a = N == 1 ? 0.5 : 0.5 - 0.25*(i-1.0)/(N-1.0)
    push!(y0, exp(t0) / (1.0 + exp(t0)) - a)
end

GC.gc()
@time sol = evolve_ode(y0, dy_dt!; parameters, precision)

#=
bytes = sizeof(sol.y) + sizeof(sol.t)
@show Base.format_bytes(bytes)
@time y, t = get_solution(sol)
=#

# BenchmarkTools.DEFAULT_PARAMETERS.seconds = 60.0
# @benchmark sol = evolve_ode(y0, dy_dt!; jacobian!, parameters)
# GC.gc()
#=
@time plot_ode(sol, method, Plots.plot)
=#

println("\ndone")
