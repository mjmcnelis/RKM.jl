
using Revise
using RKM
using LinearAlgebra
using Plots 
using UnPack
plotly()

# # TODO: need a better interface
# function dy_dt!(f, t, y)
#     A = 0.5
#     f[1] = (y[1] + A) * (1.0 - A - y[1])
# #     f .= (y .+ A) .* (1.0 .- A .- y)
#     nothing
# end

adaptive   = Fixed()
method     = Heun2()
t_span     = TimeSpan(; t0 = -10.0, tf = 10.0, dt0 = 0.001)
parameters = Parameters(; adaptive, method, t_span)

@unpack t0 = t_span
y0 = exp(t0) / (1.0 + exp(t0)) - 0.5

@time sol = evolve_ode(y0, dy_dt!; parameters)

#=
y = reduce(hcat, sol.y)'
t = sol.t 
plot(t, y; size = (900, 600), linewidth = 2,
        legend = :outertopright, legendfontsize = 12,
        ylabel = "y", yguidefontsize = 14, ytickfontsize = 12,
        xlabel = "t", xguidefontsize = 14, xtickfontsize = 12,) |> display
=#
println("\ndone")

