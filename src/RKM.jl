module RKM

import ForwardDiff: jacobian!, JacobianConfig, Dual, Chunk,
                    NANSAFE_MODE_ENABLED, DEFAULT_CHUNK_THRESHOLD
import FiniteDiff: finite_difference_jacobian!, JacobianCache
import LinearSolve: init, solve!, LinearProblem, LUFactorization,
                    AbstractFactorization#, SciMLLinearSolveAlgorithm
import LinearAlgebra: norm, dot, diagind, transpose, lu, lu!, ldiv!, mul!
import StaticArrays: SVector, SMatrix, MVector
import SparseArrays: sparse, SparseMatrixCSC
import SparseDiffTools: matrix_colors, forwarddiff_color_jacobian!, ForwardColorJacCache,
                        auto_jacvec!, DeivVecTag#, num_jacvec!
import MuladdMacro: @muladd
import FastBroadcast: @..
import UnPack: @unpack
import Base: @kwdef, rationalize, format_bytes
import DocStringExtensions: TYPEDEF, TYPEDFIELDS
import Setfield: @set!
import StatsBase: mean
# tmp for testing type stablity
import InteractiveUtils: @code_warntype

abstract type ODEMethod end

RKM_root = dirname(dirname(@__FILE__))
export RKM_root

include("timer.jl")
include("progress/default_progress.jl")
include("wrapper.jl")
include("methods/properties/embedded.jl")
include("adaptive.jl")
include("controller/pid_control.jl")
include("controller/limiter.jl")
include("controller/time_step_controller.jl")
include("implicit/root.jl")
include("implicit/jacobian.jl")
include("implicit/stage_finder.jl")
include("implicit/nansafe/nansafe_jacobian.jl")
include("implicit/nansafe/nansafe_utils.jl")
include("sensitivity/jacobian_vector.jl")
include("sensitivity/sensitivity.jl")
include("solution/solution.jl")
include("plots.jl")
include("cache.jl")
include("interpolation/interpolator_types.jl")
include("interpolation/cubic_hermite.jl")
include("interpolation/continuous_formula.jl")
include("solution/sizehint.jl")
include("solution/stats.jl")

include("methods/code_names.jl")
include("methods/properties/iteration.jl")
include("methods/properties/fesal.jl")
include("methods/properties/order.jl")
include("methods/properties/explicit_stage.jl")
# Runge-Kutta tables
include("methods/runge_kutta/runge_kutta.jl")
# include("methods/runge_kutta/debug_table.jl") # uses @test, @test_broken
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
include("interpolation/interpolation.jl")
include("sensitivity/post_sensitivity.jl")

# Adaptive methods
export Fixed, Doubling, Embedded, CentralDiff
# Time step controller
export TimeStepController, PIDControlBeta, PIDControlK
export BasicControl, PIControl, H312Control, H321PredictiveControl,
       H211bPredictiveControl
export PiecewiseLimiter, SmoothLimiter
# Implicit stage finder
export ImplicitStageFinder
# Root finder methods
export FixedPoint, Newton
# Jacobian evaluation methods
export ForwardJacobian, ForwardColorJacobian, FiniteJacobian
# Jacobian-vector evaluation methods
export NaiveJacobianVector, ForwardJacobianVector, FiniteJacobianVector
# Jacobian sparsity pattern
export nansafe_state_jacobian, nansafe_param_jacobian, test_nansafe,
       max_nan, min_nan, maximum_nan, minimum_nan
# Embedded pairs
export DefaultPair, EulerPair, SecondPair
# Numerical ODE methods
export RungeKutta, LinearMultistep
# Properties
export Iteration, Explicit, DiagonalImplicit, FullImplicit
# List of Runge-Kutta methods
export list_explicit_runge_kutta_methods, list_implicit_runge_kutta_methods
# Standard explicit Runge-Kutta
export Euler1, Heun2, Midpoint2, Ralston2, Generic2, Heun3, Ralston3, RungeKutta3,
       ShuOsher3, SpiteriRuuth3, Generic3, RungeKutta4, ThreeEightsRule4, Ralston4,
       Ketcheson4, Butcher5, Butcher6, Curtis8, Shanks8, ShanksPseudo8
# Embedded explicit Runge-Kutta
export Fehlberg12, HeunEuler21, BogackiShampine32, Fehlberg45, CashKarp54, DormandPrince54,
       BogackiShampine54, Tsitouras54, Verner56, Verner65, Fehlberg78, DormandPrince87,
       Feagin108, Feagin1210, Feagin1412
# Standard implicit Runge-Kutta
export BackwardEuler1, TrapezoidRuleBDF2, ImplicitMidpoint2,
       QinZhang2, KraaijevangerSpijker2, PareschiRusso2, PareschiRusso3, Crouzeix3,
       RadauIA3, RadauIIA3, DIRKL3, Norsett4, RaduaIA5
# Embedded implicit Runge-Kutta
export ImplicitTrapezoid21, LobattoIIIB21, LobattoIIIC21, GaussLegendre42, LobattoIIIA42,
       LobattoIIIB42, LobattoIIIC42, LobattoIIICS42, LobattoIIID42, RaduaIIA52,
       GaussLegendre64
# Linear multistep
export AdamsBashforth, AdamsMoulton, BackwardDifferentiationFormula,
       NumericalDifferentiationFormula
# Code names
export make_code_name
# Solver options
export SolverOptions
# Timer
export TimeLimit
# ODE solver
export evolve_ode, evolve_ode!
# ODE solution
export Solution, get_solution, get_stats
# Time derivatives
export get_time_derivative
# Dense output
export NoInterpolation, CubicHermite, ContinuousFormula, interpolate_solution
# Sensitivity
export NoSensitivity, DecoupledDirect, post_sensitivity_analysis, get_sensitivity
# Plots
export plot_ode
# Utilities
export debug_table

end
