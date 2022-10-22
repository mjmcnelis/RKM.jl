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
include("methods/runge_kutta/get_all_tables.jl")
        # explicit/
include("methods/runge_kutta/explicit/get_tables.jl")
            # explicit/embedded/
include("methods/runge_kutta/explicit/embedded/low_order.jl")
include("methods/runge_kutta/explicit/embedded/medium_order.jl")
include("methods/runge_kutta/explicit/embedded/high_order.jl")
include("methods/runge_kutta/explicit/embedded/very_high_order.jl")
            # explicit/fixed/
include("methods/runge_kutta/explicit/fixed/low_order.jl")
include("methods/runge_kutta/explicit/fixed/medium_order.jl")
include("methods/runge_kutta/explicit/fixed/high_order.jl")
        # implicit/
include("methods/runge_kutta/implicit/get_tables.jl")
            # implicit/fixed/
include("methods/runge_kutta/implicit/fixed/low_order.jl")
include("methods/runge_kutta/implicit/fixed/medium_order.jl")

include("parameters.jl")
include("evolve.jl")

# Adaptive methods
export Fixed, StepDoubling, Embedded, FiniteDiff
# Embedded pairs
export DefaultPair, EulerPair, SecondPair
# Numerical ODE methods
    # Properties 
export Explicit, DiagonalImplicit, FullImplicit
    # Fixed explicit Runge-Kutta
export Euler1, Heun2, Midpoint2, Ralston2, Generic2, Heun3, Ralston3, RungeKutta3, 
       ShuOsher3, SpiteriRuuth3, Generic3
export RungeKutta4, ThreeEightsRule4, Ralston4, Ketcheson4, Butcher5, Butcher6
export Curtis8, Shanks8, ShanksPseudo8
    # Embedded explicit Runge-Kutta
export Fehlberg12, HeunEuler21, BogackiShampine32
export Fehlberg45, CashKarp54, DormandPrince54, BogackiShampine54, Tsitouras54, Verner56,
       Verner65
export Fehlberg78, DormandPrince87
export Feagin108
    # Fixed implicit Runge-Kutta 
export BackwardEuler1, ImplicitMidpoint2, CrankNicolson2, QinZhang2, KraaijevangerSpijker2,
       PareschiRusso2, LobattoIIIB2, LobattoIIIC2, PareschiRusso3, Crouzeix3, 
       RadauIA3, RadauIIA3, DIRKL3
export Norsett4, LobattoIIIA4, LobattoIIIB4, LobattoIIIC4
# Get Runge-Kutta tables
export get_all_runge_kutta_tables, get_runge_kutta_explicit_tables, 
       get_runge_kutta_full_implicit_tables, get_runge_kutta_diagonal_implicit_tables
# Code names
export get_code_name
# ODE solution
export Solution
# ODE options
export Parameters
# ODE solver
export evolve_ode
# Utilities 
export debug_table, debug_iteration

end