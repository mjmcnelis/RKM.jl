module RKM

using DataFrames
using UnPack
import Base.@kwdef

abstract type DiffEqMethod end

include("utils.jl")
include("solution.jl")
include("embedded.jl")
include("adaptive.jl")

include("methods/iteration.jl")
include("methods/runge_kutta/runge_kutta.jl")
include("methods/runge_kutta/explicit/fixed/low_order.jl")
include("methods/runge_kutta/explicit/fixed/medium_order.jl")


include("parameters.jl")
include("evolve.jl")

# Adaptive methods
export Fixed, StepDoubling, Embedded, FiniteDiff
# Embedded pairs
export DefaultPair, EulerPair, SecondPair
# Numerical methods (explicit Runge-Kutta)
export Euler1, Heun2, Midpoint2, Ralston2
export Heun3, Ralston3, RungeKutta3, ShuOsher3, SpiteriRuuth3
export RungeKutta4, ThreeEightsRule4, Ralston4, Ketcheson4
# ODE solution
export Solution
# ODE options
export Parameters
# ODE solver
export evolve_ode

end