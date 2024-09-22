using Revise, RKM
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/logistic/equations.jl") : nothing
include("$RKM_root/validation/ode/logistic/parameters.jl")

options = SolverOptions(options)

t0 = -5.0
tf = 5.0
dt0 = 0.1

N = 1000
p = [0.5 - 0.25*(i-1.0)/(N-1.0+eps(1.0)) for i in 1:N]
y0 = Float64[]
for i = eachindex(p)
    push!(y0, exp(t0) / (1.0 + exp(t0)) - p[i])
end

@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options; model_parameters = p)

@time yp = post_sensitivity_analysis(sol, options, dy_dt!, p)
# @show Base.format_bytes(sizeof(yp))

ny = sol.dimensions[1]
nt = length(sol.t)
np = length(p)

# time slice: dy_1/dp1 dy2/dp1 dy1/dp2 dy2/dp2
yp = reshape(yp, ny*np, nt) |> transpose

# get time slice (transpose?)
# @time a = reshape(view(yp, 2, :), ny, np)

# get parameter slice
# @time b = view(yp, :, 1:ny)

# plot(sol.t, yp, ylims = (0.0, 50.0)) |> display

# get_stats(sol)

GC.gc()
println("\ndone")
