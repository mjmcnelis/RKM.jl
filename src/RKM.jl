module RKM

import ProgressMeter: Progress, next!
import LinearAlgebra: norm, tril, diag
import StaticArrays: SVector, SMatrix, MVector, MMatrix, @MVector, @MMatrix
import MuladdMacro: @muladd
import FastBroadcast: @..
import UnPack: @unpack
import Base: @kwdef
import Test: @test, @test_broken
import DocStringExtensions: TYPEDEF, TYPEDFIELDS

abstract type ODEMethod end

const VectorMVector = Union{Vector{T}, MVector{D,T}} where {D, T <: AbstractFloat} 
const MatrixMMatrix = Union{Matrix{T}, MMatrix{D,S,T,DS}} where {D, S, DS, 
                                                                 T <: AbstractFloat}

RKM_root = dirname(dirname(@__FILE__))
export RKM_root

include("time.jl")
include("solution.jl")
include("embedded.jl")
include("adaptive.jl")
include("plots.jl")
include("utils.jl")

include("methods/code_names.jl")
include("methods/properties.jl")
include("methods/utils.jl")
# Runge-Kutta tables
include("methods/runge_kutta/runge_kutta.jl")
include("methods/runge_kutta/debug_table.jl")
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
include("updates/function_evaluations.jl")
include("updates/step_rejection_rate.jl")
include("updates/jacobian.jl")
include("updates/runge_kutta/explicit/update.jl")
include("updates/runge_kutta/implicit/update.jl")

include("parameters.jl")
include("evolve.jl")

# Adaptive methods
export Fixed, Doubling, Embedded, FiniteDiff
# Embedded pairs
export DefaultPair, EulerPair, SecondPair
# Numerical ODE methods
export RungeKutta
    # Properties
export Iteration, Explicit, DiagonalImplicit, FullImplicit, 
       FirstSameAsLast, FSAL, NotFSAL
    # Fixed explicit Runge-Kutta
export Euler1, Heun2, Midpoint2, Ralston2, Generic2, Heun3, Ralston3, RungeKutta3,
       ShuOsher3, SpiteriRuuth3, Generic3, RungeKutta4, ThreeEightsRule4, Ralston4, 
       Ketcheson4, Butcher5, Butcher6, Curtis8, Shanks8, ShanksPseudo8
    # Embedded explicit Runge-Kutta
export Fehlberg12, HeunEuler21, BogackiShampine32, Fehlberg45, CashKarp54, DormandPrince54, 
       BogackiShampine54, Tsitouras54, Verner56, Verner65, Fehlberg78, DormandPrince87,
       Feagin108, Feagin1210, Feagin1412
    # Fixed implicit Runge-Kutta
export BackwardEuler1, ImplicitMidpoint2, QinZhang2, KraaijevangerSpijker2, PareschiRusso2,
       PareschiRusso3, Crouzeix3, RadauIA3, RadauIIA3, DIRKL3, Norsett4, RaduaIA5
    # Embedded implicit Runge-Kutta
export CrankNicolson21, LobattoIIIB21, LobattoIIIC21, GaussLegendre42, LobattoIIIA42, 
       LobattoIIIB42, LobattoIIIC42, LobattoIIICS42, LobattoIIID42, RaduaIIA52, GaussLegendre64
# Code names
export make_code_name
# ODE solution
export get_solution
# Parameters
export Parameters
# Time
export TimeRange, TimeLimit
# ODE solver
export evolve_ode
# Plots
export plot_ode
# Utilities
export debug_table

end
