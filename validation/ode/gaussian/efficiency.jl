using Revise, RKM, LinearAlgebra, OrdinaryDiffEq
import DoubleFloats: Double64
import StaticArrays: SA
import BenchmarkTools: @benchmark, @btime, mean
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/gaussian/equations.jl") : nothing

precision = Float64
# precision = Double64

t0 = -5.0 |> precision
y0 = [exp(-t0^2/2.0) + 1.0]
dt0 = 1e-4
tf = 5.0

safety = 0.9
low = 0.1
high = 5.0

fe_ord = Float64[]
fe_rkm = Float64[]
err_ord = Float64[]
err_rkm = Float64[]

epsilon_vect = 10.0.^(-(5:1:15)) # Float64
# epsilon_vect = 10.0.^(-(5:1:26)) # Double64

for epsilon in epsilon_vect
    @show epsilon
    # OrdinaryDiffEq
    prob = ODEProblem(f_ord, y0, (t0, tf))
    sol = solve(prob, DP5(), dt = dt0, reltol = epsilon, abstol = 0.0,
                qmin = 1/high, qmax = 1/low, gamma = safety, controller = IController())

    append!(fe_ord, sol.destats.nf)

    y, t = vcat(sol.u...), sol.t
    y_ex = zeros(BigFloat, size(y)...)
    err = zeros(BigFloat, length(t))
    for i in eachindex(t)
        i == 1 ? continue : nothing
        y_ex[i,:] = y_exact(t[i])
        err[i] = norm(y[i,:] .- y_ex[i,:], 2)
    end
    append!(err_ord, mean(err))

    # RKM
    ps = SolverOptions(; adaptive = Embedded(; epsilon, low, high, safety),
                      method = DormandPrince5(), t_range = TimeRange(; t0, tf, dt0))
    sol = evolve_ode(y0, dy_dt!; options = ps, show_progress = false, precision)

    append!(fe_rkm, sol.FE)

    t, y = get_solution(sol)
    y_ex = zeros(BigFloat, size(y)...)
    err = zeros(BigFloat, length(t))
    for i in eachindex(t)
        i == 1 ? continue : nothing
        y_ex[i,:] = y_exact(t[i])
        err[i] = norm(y[i,:] .- y_ex[i,:], 2)
    end
    append!(err_rkm, mean(err))
end

plot(; legend = :outertopright, legendtitlefontsize = 12, legendfontsize = 12,
       size = (900, 600), xlabel = "Function evaluations", ylabel = "Mean error",
       yguidefontsize = 14, ytickfontsize = 12, yaxis = :log, ylims = (1e-30, 1e0),
       xguidefontsize = 14, xtickfontsize = 12, xaxis = :log, xlims = (1e2, 1e7),
    )
plot!(fe_ord, err_ord, linewidth = 2, label = "DP5 OrdinaryDiffEq")
plot!(fe_rkm, err_rkm, linewidth = 2, label = "DP54 RKM")
