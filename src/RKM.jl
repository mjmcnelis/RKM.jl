module RKM

import Dates: now, Minute, DateTime
import LinearAlgebra: norm, tril, diag
import Statistics: mean
import DoubleFloats: Double64
import StaticArrays: SVector, SMatrix, MVector
import DataStructures: OrderedDict
import MuladdMacro: @muladd
# note: fast broadcasting doesn't work with BigFloat MVector (t, dt)
import FastBroadcast: @..
import UnPack: @unpack
import Setfield: @set
import Base: @kwdef

abstract type ODEMethod end

RKM_root = dirname(dirname(@__FILE__))

include("time.jl")
include("solution.jl")
include("embedded.jl")
include("adaptive.jl")
include("plots.jl")

include("methods/code_names.jl")
include("methods/properties.jl")
include("methods/utils.jl")
# Runge-Kutta tables
include("methods/runge_kutta/runge_kutta.jl")
include("methods/runge_kutta/explicit/fixed/low_order.jl")
include("methods/runge_kutta/explicit/fixed/medium_order.jl")
include("methods/runge_kutta/explicit/fixed/high_order.jl")
include("methods/runge_kutta/explicit/embedded/low_order.jl")
include("methods/runge_kutta/explicit/embedded/medium_order.jl")
include("methods/runge_kutta/explicit/embedded/high_order.jl")
include("methods/runge_kutta/explicit/embedded/very_high_order.jl")
include("methods/runge_kutta/implicit/fixed/low_order.jl")
include("methods/runge_kutta/implicit/fixed/medium_order.jl")
include("methods/runge_kutta/implicit/embedded/low_order.jl")
include("methods/runge_kutta/implicit/embedded/medium_order.jl")
# Runge-Kutta updates
include("updates/evaluations.jl")
include("updates/jacobian.jl")
include("updates/runge_kutta/explicit/update.jl")
include("updates/runge_kutta/implicit/update.jl")

include("parameters.jl")
include("evolve.jl")
include("efficiency.jl")

# Adaptive methods
export Fixed, Doubling, Embedded, FiniteDiff
# Embedded pairs
export DefaultPair, EulerPair, SecondPair
# Numerical ODE methods
export RungeKutta
    # Properties
export Iteration, Explicit, DiagonalImplicit, FullImplicit
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
export BackwardEuler1, ImplicitMidpoint2, QinZhang2, KraaijevangerSpijker2, PareschiRusso2,
       PareschiRusso3, Crouzeix3, RadauIA3, RadauIIA3, DIRKL3
export Norsett4, RaduaIA5
    # Embedded implicit Runge-Kutta
export CrankNicolson21, LobattoIIIB21, LobattoIIIC21
export GaussLegendre42, LobattoIIIA42, LobattoIIIB42, LobattoIIIC42, LobattoIIICS42,
       LobattoIIID42, RaduaIIA52, GaussLegendre64
# Get Runge-Kutta tables
export get_all_runge_kutta_tables, get_runge_kutta_explicit_tables,
       get_runge_kutta_full_implicit_tables, get_runge_kutta_diagonal_implicit_tables
# Code names
export make_code_name
# ODE solution
export Solution, get_solution
# Parameters
export Parameters
# Time
export TimeRange, TimeLimit
# ODE solver
export evolve_ode
# Efficiency 
export efficiency_curve
# Utilities
export debug_table, debug_iteration, RKM_root
# Plots
export plot_ode

end
