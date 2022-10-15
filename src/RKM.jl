module RKM

# using DataFrames      # may use to output RK tables
using LinearAlgebra
using Test
using UnPack
import Base.@kwdef

abstract type ODEMethod end

include("utils.jl")
include("solution.jl")
include("embedded.jl")
include("adaptive.jl")

# methods/
include("methods/properties.jl")
    # runge_kutta/
include("methods/runge_kutta/runge_kutta.jl")
include("methods/runge_kutta/utils.jl")
       # explicit/embedded/
include("methods/runge_kutta/explicit/embedded/low_order.jl")
include("methods/runge_kutta/explicit/embedded/medium_order.jl")
include("methods/runge_kutta/explicit/embedded/high_order.jl")
include("methods/runge_kutta/explicit/embedded/very_high_order.jl")
        # explicit/fixed/
include("methods/runge_kutta/explicit/fixed/low_order.jl")
include("methods/runge_kutta/explicit/fixed/medium_order.jl")
include("methods/runge_kutta/explicit/fixed/high_order.jl")

include("parameters.jl")
include("evolve.jl")

# Adaptive methods
export Fixed, StepDoubling, Embedded, FiniteDiff
# Embedded pairs
export DefaultPair, EulerPair, SecondPair
# Numerical methods (fixed explicit Runge-Kutta)
export Euler1, Heun2, Midpoint2, Ralston2, Generic2
export Heun3, Ralston3, RungeKutta3, ShuOsher3, SpiteriRuuth3, Generic3
export RungeKutta4, ThreeEightsRule4, Ralston4, Ketcheson4
export Butcher5, Butcher6
export Curtis8, Shanks8, ShanksPseudo8
# Numerical methods (embedded explicit Runge-Kutta)
export Fehlberg12, HeunEuler21, BogackiShampine32
export Fehlberg45, CashKarp54, DormandPrince54, BogackiShampine54
export Tsitouras54, Verner56, Verner65
export Fehlberg78, DormandPrince87
export Feagin108
# Code names
export get_code_name
# ODE solution
export Solution
# ODE options
export Parameters
# ODE solver
export evolve_ode
# Utilities 
export debug_table

end