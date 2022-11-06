
using Revise
using RKM
using LinearAlgebra
using Plots 
using UnPack
plotly()
# reduces allocations but can't redefine it
try
    dy_dt!
catch err
    isa(err, UndefVarError) ? include(joinpath(RKM_root, "equations.jl")) : nothing
end

# adaptive   = Fixed()
adaptive   = Doubling()
method     = Heun2()
t_span     = TimeSpan(; t0 = -10.0, tf = 10.0, dt0 = 0.001)

parameters = Parameters(; adaptive, method, t_span)

@unpack t0 = t_span
y0 = exp(t0) / (1.0 + exp(t0)) - 0.5

# TEMP so can see how solver works with vectors
y02 = exp(t0) / (1.0 + exp(t0)) - 0.25
y0 = [y0, y02]

@time sol = evolve_ode(y0, dy_dt!; parameters)

function plot_y(sol)
    y = reduce(hcat, sol.y)'
    t = sol.t 
    # TODO: may have Solution outer constructor to compute this
    dt = sol.t[2:end] .- sol.t[1:end-1]
    plot(t, y; size = (900, 600), linewidth = 2,
         legend = :outertopright, legendfontsize = 12,
         ylabel = "y", yguidefontsize = 14, ytickfontsize = 12,
         xlabel = "t", xguidefontsize = 14, xtickfontsize = 12,) |> display
    plot(t[1:end-1], dt; size = (900, 600), linewidth = 2,
         label = "", legend = :outertopright, legendfontsize = 12,
         ylabel = "Δt", yguidefontsize = 14, ytickfontsize = 12,
         xlabel = "t", xguidefontsize = 14, xtickfontsize = 12,) |> display
end
#=
plot_y(sol)
=#

println("\ndone")

