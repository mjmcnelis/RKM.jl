using Revise, RKM, BenchmarkTools
import DoubleFloats: Double64
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/exponential/equations.jl") : nothing

# p = [1.0]           # exponential growth
p = [-1.0]          # exponential decay

t0 = 0.0 |> BigFloat
y0 = exp(p[1] * t0)

tf = 10.0
dt0 = 1e-4

options = SolverOptions(method = BackwardEuler1(),
                        adaptive = Fixed(),
                        sensitivity_method = DecoupledDirect()
                       )

@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p)

t, y = get_solution(sol)
t, S = get_sensitivity(sol)

get_stats(sol)

plot_ode(sol, options.method, Plots.plot)#; logy = true, show_time_step = true)
plot!(sol.t, y_exact.(sol.t; p), label = "exact") |> display

# exact 1st order sensitivity (S)
f1(t) = t*exp(p[1]*t)
plot(f1, t0, tf; color = :red, linewidth = 2)

# numerical solution
plot!(t, S; color = :black, line = :dash, linewidth = 2) |> display

# exact 2nd order sensitivity (R)
f2(t) = t^2*exp(p[1]*t)
plot(f2, t0, tf; color = :red, linewidth = 2)#, yaxis = :log)

# finite differences
z(t,p) = exp(p*t)
dp = 1e-5
f2_FD(t) = (z(t,p[1]+dp) - 2*z(t,p[1]) + z(t,p[1]-dp))/dp^2
plot!(f2_FD, t0, tf; color = :black, linewidth = 2, line = :dash)

S_func(y; p) = y*log(y)/p

yP = y .+ dp*S
pP = p[1] + dp

yM = y .- dp*S
pM = p[1] - dp

plot!(t, (S_func.(yP; p = pP) .- S_func.(y; p = p[1]))./(dp);
      color = :blue, line = :dot, linewidth = 2)

GC.gc()

solP = deepcopy(sol)
solP.y .= y .+ vec(dp*S)
SP = post_sensitivity_analysis(solP, options, dy_dt!, p .+ dp)

solM = deepcopy(sol)
solM.y .= y .- vec(dp*S)
SM = post_sensitivity_analysis(solM, options, dy_dt!, p .- dp)

# why sqrt(2) and not 2?
yPP = yP .+ dp*SP
yMM = yM .- dp*SM
R = (yPP .- 2.0.*y .+ yMM) ./ (sqrt(2)*dp)^2
# R = (yPP .- 2.0.*yP .+ y) ./ (dp)^2
# R = (y .- 2.0.*yM .+ yMM) ./ (dp)^2

# TODO: why does forward/backward difference do worse with smaller dp
# resolved: needed to use BackwardEuler1 for original sensitivity S
# R = (SP .- SM) ./ (2*dp)
# R = (SP .- S) ./ dp
# R = (S .- SM) ./ dp
plot!(t, R; color = :green, linewidth = 2, line = :dot) |> display
println("\ndone")

# plot(t, S_func.(yP; p = pP))
# plot!(t, SP)

# plot(t, S_func.(yM; p = pM))
# plot!(t, SM)