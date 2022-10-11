module RKM

using UnPack
import Base.@kwdef

abstract type DiffEqMethod end

include("solution.jl")
include("embedded.jl")
include("adaptive.jl")

include("methods/iteration.jl")
include("methods/runge_kutta/runge_kutta.jl")
include("methods/runge_kutta/explicit/fixed/low_order.jl")


include("parameters.jl")
include("evolve.jl")

# Adaptive methods
export Fixed, StepDoubling, Embedded, FiniteDiff
# Embedded pairs
export DefaultPair, EulerPair, SecondPair
# Numerical methods
export Euler1, Heun2
# ODE solution
export Solution
# ODE options
export Parameters
# ODE solver
export evolve_ode

end