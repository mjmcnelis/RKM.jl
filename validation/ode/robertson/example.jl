using Revise, RKM, BenchmarkTools
import DoubleFloats: Double64
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/rober/equations.jl") : nothing

# adaptive = Fixed()   
adaptive = Doubling(; epsilon = 1e-6)  

method = TrapezoidRuleBDF2()  
# method = BackwardEuler1()

controller = PIDControllerK(; kI = 0.3, kP = 0.4)

# stage_finder = ImplicitStageFinder()
stage_finder = ImplicitStageFinder(; jacobian_method = ForwardJacobian())

t_range = TimeRange(; t0 = 0.01, tf = 1.0e4, dt0 = 0.01)

parameters = Parameters(; adaptive, method, t_range, controller, stage_finder)

y0 = [1.0, 0.0, 0.0]

@time sol = evolve_ode(y0, dy_dt!; parameters)

get_stats(sol)
# plot_ode(sol, method, Plots.plot; logx = true)

GC.gc()
println("\ndone")
