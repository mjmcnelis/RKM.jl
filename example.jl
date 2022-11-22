
using Revise
using RKM
using LinearAlgebra
using Plots 
using UnPack
plotly()
# reduces allocations but can't redefine it
try
    dy_dt!
    jacobian!
catch err
    isa(err, UndefVarError) ? include("$RKM_root/equations.jl") : nothing
end

# TODO: try to see if I can rescale epsilon? 

adaptive = Fixed()
# adaptive = Doubling()
# adaptive = Embedded()

# method = BackwardEuler1()
# method = Euler1()
method = Heun2()
# method = HeunEuler21()
# method = Fehlberg45()

t_span     = TimeSpan(; t0 = -10.0, tf = 10.0, dt0 = 0.001)

# do asserts between adaptive, method in parameters outer-constructor
parameters = Parameters(; adaptive, method, t_span)

@unpack t0 = t_span
y0 = exp(t0) / (1.0 + exp(t0)) - 0.5

# TEMP so can see how solver works with vectors
y02 = exp(t0) / (1.0 + exp(t0)) - 0.25
y0 = [y0, y02]

@time sol = evolve_ode(y0, dy_dt!; jacobian!, 
                                   parameters)

#=
plot_ode(sol, method, Plots.plot)
=#

println("\ndone")

