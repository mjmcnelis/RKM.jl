
using Revise
using RKM
using LinearAlgebra
using Plots 
using UnPack
plotly()
# reduces allocations but can't redefine it
try
    dy_dt!
catch err
    isa(err, UndefVarError) ? include(joinpath(RKM_root, "equations.jl")) : nothing
end

# TODO: try to see if I can rescale epsilon? 

# adaptive   = Fixed()
adaptive   = Doubling()
method     = Heun2()
t_span     = TimeSpan(; t0 = -10.0, tf = 10.0, dt0 = 0.001)

parameters = Parameters(; adaptive, method, t_span)

@unpack t0 = t_span
y0 = exp(t0) / (1.0 + exp(t0)) - 0.5

# TEMP so can see how solver works with vectors
y02 = exp(t0) / (1.0 + exp(t0)) - 0.25
y0 = [y0, y02]

@time sol = evolve_ode(y0, dy_dt!; parameters)

#=
plot_ode(sol, method, Plots.plot)
=#

println("\ndone")

