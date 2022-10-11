
using Revise
using RKM

adaptive = Embedded()
method   = Heun2()

parameters = Parameters(; adaptive, method)

@show method


q()

sol = evolve_ode(; parameters)