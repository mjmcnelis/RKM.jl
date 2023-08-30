using Revise, RKM, BenchmarkTools
import DoubleFloats: Double64
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/logistic/equations.jl") : nothing

precision = Float64
# precision = Double64
# precision = BigFloat        # 31.60 M allocations (fixed time step, no progress meter)

# adaptive = Fixed()
# method = RungeKutta4()
# method = BackwardEuler1()     # 200.31 k allocations: 19.856 MiB

# TODO: if stage 2 and 3 have "same jacobian", then can reuse it?
#       not entirely sure what it implies
# method = TrapezoidRuleBDF2()  # 400.41 k allocations: 35.122 MiB (fixed time step)
# method = CrankNicolson21()    # 200.34 k allocations: 19.858 MiB (fixed time step)
# method = Heun2()

# adaptive = Doubling()
adaptive = Embedded()#; epsilon = 1e-7)
method = HeunEuler21()
# method = Fehlberg45()
# adaptive = FiniteDiff()
# method = Heun2()

# pid = BasicControl()
# pid = PIControl()
pid = H312Control()
# pid = H321PredictiveControl()
# pid = H211bPredictiveControl()

controller = TimeStepController(; pid)

stage_finder = ImplicitStageFinder()
# stage_finder = ImplicitStageFinder(; root_method = FixedPoint())
# stage_finder = ImplicitStageFinder(; jacobian_method = ForwardJacobian())

t_range = TimeRange(; t0 = -10, tf = 10, dt0 = 1e-4)

# do asserts between adaptive, method in parameters outer-constructor
parameters = Parameters(; adaptive, method, t_range, controller, stage_finder)

t0 = t_range.t0

N = 2
y0 = Float64[]
for i = 1:N
    local a = N == 1 ? 0.5 : 0.5 - 0.25*(i-1.0)/(N-1.0)
    push!(y0, exp(t0) / (1.0 + exp(t0)) - a)
end

# show_progress = true
show_progress = false

static_array = false
# static_array = true

@time sol = evolve_ode(y0, dy_dt!; parameters, precision, show_progress, static_array)

# in-place version
# sol = Solution()
# @time evolve_ode!(sol, y0, dy_dt!; parameters, show_progress, static_array)

get_stats(sol)
# plot_ode(sol, method, Plots.plot)

GC.gc()
println("\ndone")
