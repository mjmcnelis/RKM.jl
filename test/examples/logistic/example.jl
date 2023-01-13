using Revise, RKM, BenchmarkTools
import DoubleFloats: Double64
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/test/examples/logistic/equations.jl") : nothing 

precision = Float64
# precision = Double64
# precision = BigFloat        # 31.60 M allocations (fixed time step)

adaptive = Fixed()         
method = RungeKutta4()
# adaptive = Doubling()        
# method = Heun2()
# adaptive = Embedded()       
# method = Fehlberg45()

t_range = TimeRange(; t0 = -10, tf = 10, dt0 = 1e-4)

# do asserts between adaptive, method in parameters outer-constructor
parameters = Parameters(; adaptive, method, t_range)

t0 = t_range.t0

N = 2
y0 = Float64[]
for i = 1:N
    local a = N == 1 ? 0.5 : 0.5 - 0.25*(i-1.0)/(N-1.0)
    push!(y0, exp(t0) / (1.0 + exp(t0)) - a)
end

show_progress = true
# show_progress = false

@time sol = evolve_ode(y0, dy_dt!; parameters, precision, show_progress)
# @btime sol = evolve_ode(y0, dy_dt!; parameters, precision, show_progress)

# @show Base.format_bytes(sizeof(sol.y) + sizeof(sol.t))
# plot_ode(sol, method, Plots.plot)

GC.gc()
println("\ndone")
