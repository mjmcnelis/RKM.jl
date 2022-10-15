
using Revise
using RKM

adaptive = Embedded()
method   = Heun2()

Euler1(), Heun2(), Midpoint2(), Ralston2()
export Heun3(), Ralston3(), RungeKutta3(), ShuOsher3(), SpiteriRuuth3()
export RungeKutta4(), ThreeEightsRule4(), Ralston4(), Ketcheson4()

parameters = Parameters(; adaptive, method)

@show method


q()

sol = evolve_ode(; parameters)