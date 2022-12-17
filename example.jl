using Revise
using RKM
using LinearAlgebra
using StaticArrays
using Plots 
using UnPack
plotly()
try
    dy_dt!
    jacobian!
catch err
    isa(err, UndefVarError) ? include("$RKM_root/equations.jl") : nothing
end
# TODO: try to see if I can rescale epsilon? 

adaptive = Fixed()             # 200k allocs (dt = 1e-4)
# adaptive = Doubling()          # 281 allocs (Heun2)
# adaptive = Embedded()          # 60 allocs (Fehlberg45)

method = RungeKutta4()
# method = Heun2()
# method = Fehlberg45()
# method = BackwardEuler1()

t_span = TimeSpan(; t0 = -10.0, tf = 10.0, dt0 = 1e-4)
timer = TimeLimit()
 
data_format = TimeSlice()
# data_format = SpaceSlice()

# do asserts between adaptive, method in parameters outer-constructor
parameters = Parameters(; adaptive, method, t_span, timer, data_format)

@unpack t0 = t_span

N = 2
y0 = Float64[]
for i = 1:N
    global a = N == 1 ? 0.5 : 0.5 - 0.25*(i-1.0)/(N-1.0)
    push!(y0, exp(t0) / (1.0 + exp(t0)) - a)
end

GC.gc()
@time sol = evolve_ode(y0, dy_dt!; jacobian!, parameters)
GC.gc()
#=
plot_ode(sol, method, Plots.plot)
=#

println("\ndone")
