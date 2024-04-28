module RKM

import SciMLBase: init, solve!
import ForwardDiff: jacobian!, JacobianConfig, DEFAULT_CHUNK_THRESHOLD
import FiniteDiff: finite_difference_jacobian!, JacobianCache
import LinearSolve: LinearProblem, LinearCache, set_A, set_b, do_factorization,
                    set_cacheval, _ldiv!, KLUFactorization, LUFactorization,
                    AbstractFactorization
import LinearAlgebra
import LinearAlgebra: norm, tril, diag, diagind, lu, lu!, eigvals, transpose
import StaticArrays: SVector, SMatrix, MVector, MMatrix, @MVector, @MMatrix
import SparseArrays: SparseMatrixCSC, sparse
import MuladdMacro: @muladd
import FastBroadcast: @..
import DoubleFloats: DoubleFloat, IEEEFloat
import UnPack: @unpack
import Base: @kwdef, rationalize, format_bytes
import ProgressMeter: Progress, next!
import Test: @test, @test_broken
import DocStringExtensions: TYPEDEF, TYPEDFIELDS
import Setfield: @set!
import InteractiveUtils: @code_warntype

abstract type ODEMethod end

RKM_root = dirname(dirname(@__FILE__))
export RKM_root

include("timer.jl")
include("progress.jl")
include("wrapper.jl")
include("embedded.jl")
include("adaptive.jl")
include("tmp/double_float.jl")
include("controller/pid_control.jl")
include("controller/limiter.jl")
include("controller/time_step_controller.jl")
include("stage_finder.jl")
include("solution.jl")
include("tmp/linear_solver.jl")
include("plots.jl")
include("cache.jl")
include("interpolation.jl")

include("methods/code_names.jl")
include("methods/properties.jl")
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
# Multistep tables
include("methods/multistep/linear_multistep.jl")
include("methods/multistep/adams/adams_bashforth.jl")
include("methods/multistep/adams/adams_moulton.jl")
include("methods/multistep/differentiation_formula/backward_differentiation_formula.jl")
include("methods/multistep/differentiation_formula/numerical_differentiation_formula.jl")
# Reconstruct ODE method
include("methods/utils.jl")
# Runge-Kutta updates
include("updates/runge_kutta/runge_kutta_step.jl")
include("updates/runge_kutta/fixed_step.jl")
include("updates/runge_kutta/double_step.jl")
include("updates/runge_kutta/embedded_step.jl")
include("updates/runge_kutta/central_step.jl")
# Multistep updates
include("updates/multistep/adams/adams_step.jl")
include("updates/multistep/adams/fixed_step.jl")

include("options.jl")
include("evolve.jl")

# Adaptive methods
export Fixed, Doubling, Embedded, CentralDiff
# Time step controller
export TimeStepController, PIDControlBeta, PIDControlK
export BasicControl, PIControl, H312Control, H321PredictiveControl,
       H211bPredictiveControl
export PiecewiseLimiter, SmoothLimiter
# Implicit stage finder
export ImplicitStageFinder, FixedPoint, Newton, ForwardJacobian, FiniteJacobian
# Embedded pairs
export DefaultPair, EulerPair, SecondPair
# Interpolation methods
export NoInterpolator, HermiteInterpolator

# TODO: probably don't need to export these
#---------------------------------------
# Numerical ODE methods
export RungeKutta, LinearMultistep
    # Properties
export Iteration, Explicit, DiagonalImplicit, FullImplicit,
       FirstSameAsLast, FSAL, NotFSAL
#---------------------------------------
    # Standard explicit Runge-Kutta
export list_explicit_runge_kutta_methods, list_implicit_runge_kutta_methods
export Euler1, Heun2, Midpoint2, Ralston2, Generic2, Heun3, Ralston3, RungeKutta3,
       ShuOsher3, SpiteriRuuth3, Generic3, RungeKutta4, ThreeEightsRule4, Ralston4,
       Ketcheson4, Butcher5, Butcher6, Curtis8, Shanks8, ShanksPseudo8
    # Embedded explicit Runge-Kutta
export Fehlberg12, HeunEuler21, BogackiShampine32, Fehlberg45, CashKarp54, DormandPrince54,
       BogackiShampine54, Tsitouras54, Verner56, Verner65, Fehlberg78, DormandPrince87,
       Feagin108, Feagin1210, Feagin1412
    # Standard implicit Runge-Kutta
export BackwardEuler1, ImplicitMidpoint2, QinZhang2, KraaijevangerSpijker2, PareschiRusso2,
       PareschiRusso3, Crouzeix3, RadauIA3, RadauIIA3, DIRKL3, Norsett4, RaduaIA5
# not sure which category it belongs to yet
export TrapezoidRuleBDF2
    # Embedded implicit Runge-Kutta
export CrankNicolson21, LobattoIIIB21, LobattoIIIC21, GaussLegendre42, LobattoIIIA42,
       LobattoIIIB42, LobattoIIIC42, LobattoIIICS42, LobattoIIID42, RaduaIIA52, GaussLegendre64
    # Linear multistep
export AdamsBashforth, AdamsMoulton, BackwardDifferentiationFormula,
       NumericalDifferentiationFormula

# Code names
export make_code_name
# ODE solution
export Solution, get_solution, clear_solution!, get_stats
# Solver options
export SolverOptions
# Time
export TimeRange, TimeLimit
# ODE solver
export evolve_ode, evolve_ode!
# Plots
export plot_ode
# Utilities
export debug_table

end
