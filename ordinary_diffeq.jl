using Revise
using OrdinaryDiffEq
using StaticArrays
using Plots
plotly()
# function fp(f, y, p, t)
#     f[1] = (y[1] + 0.5) * (0.5 - y[1])
#     f[2] = (y[2] + 0.25) * (0.75 - y[2])
#     nothing
# end
# function fp(y, p, t)
#      SA[(y[1] + 0.5) * (0.5 - y[1]),
#         (y[2] + 0.25) * (0.75 - y[2])]
# end
t0 = -10.0
# u0 = [exp(t0) / (1.0 + exp(t0)) - 0.5,
#       exp(t0) / (1.0 + exp(t0)) - 0.25]
u0 = SA[exp(t0) / (1.0 + exp(t0)) - 0.5;
        exp(t0) / (1.0 + exp(t0)) - 0.25]

tspan = (t0, 10.0)
prob = ODEProblem(fp, u0, tspan)

GC.gc()
integ = init(prob, Euler(), dt = 5e-6)
@time solve!(integ)
sol = integ.sol

# @show sol.destats
println("\ndone")

#=
plot(sol,linewidth=5,title="Solution to the linear ODE with a thick line",
     xaxis="Time (t)",yaxis="u(t) (in μm)",label="My Thick Line!") # legend=false
=#