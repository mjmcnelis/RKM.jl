using Revise, RKM, BenchmarkTools
import DoubleFloats: Double64
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/sine/equations.jl") : nothing 

const ω = 100.0       # frequency
const R = 20          # cycles

precision = Float64
adaptive = Fixed()         
method = RungeKutta4()
t_range = TimeRange(; t0 = 0, tf = pi/2, dt0 = 1e-5)
parameters = Parameters(; adaptive, method, t_range)

t0 = t_range.t0 |> precision
y0 = [sin(t0), ω*cos(ω*t0)]

@time sol = evolve_ode(y0, dy_dt!; parameters, precision)
# @btime sol = evolve_ode(y0, dy_dt!; parameters, precision)
# @show Base.format_bytes(sizeof(sol.y) + sizeof(sol.t))
# plot_ode(sol, method, Plots.plot)

GC.gc()
println("\ndone")
