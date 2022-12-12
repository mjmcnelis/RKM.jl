# note: used to benchmark performance against OrdinaryDiffEq
using Revise
using OrdinaryDiffEq
import RKM: RKM_root
using StaticArrays
using Plots
plotly()
try
    fp
catch err
    isa(err, UndefVarError) ? include("$RKM_root/equations.jl") : nothing
end

t0 = -10.0
# u0 = [exp(t0) / (1.0 + exp(t0)) - 0.5,
#       exp(t0) / (1.0 + exp(t0)) - 0.25]
u0 = SA[exp(t0) / (1.0 + exp(t0)) - 0.5;
        exp(t0) / (1.0 + exp(t0)) - 0.25]

tspan = (t0, 10.0)
prob = ODEProblem(fp, u0, tspan)

GC.gc()
integ = init(prob, RK4(), dt = 1e-4, adaptive = false, saveat = 1000)
@time solve!(integ)
sol = integ.sol
GC.gc()

# @show sol.destats
println("\ndone")

#=
plot(sol,linewidth=5,title="Solution to the linear ODE with a thick line",
     xaxis="Time (t)",yaxis="u(t) (in μm)",label="My Thick Line!") # legend=false
=#