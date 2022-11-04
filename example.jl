
using Revise
using RKM
using LinearAlgebra
using Plots 
using UnPack
plotly()
# TODO: wrap into function
# function load_equations(file_path)
#     try
#         dy_dt!    # but this is a global...
#     catch err
#         isa(err, UndefVarError) ? include(file_path) : nothing
#         return dy_dt!
#     end
# end
# file_path = joinpath(RKM_root, "equations.jl")
# dy_dt! = load_equations(file_path)
try
    dy_dt!
catch err
    isa(err, UndefVarError) ? include(joinpath(RKM_root, "equations.jl")) : nothing
end

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

