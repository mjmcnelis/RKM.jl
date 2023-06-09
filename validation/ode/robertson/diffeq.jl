using Revise, OrdinaryDiffEq, StaticArrays, BenchmarkTools
using DoubleFloats: Double64
import RKM: RKM_root
using Plots; plotly()
!(@isdefined f_ord) ? include("$RKM_root/validation/ode/logistic/equations.jl") : nothing 

y0 = [1.0, 0.0, 0.0]
alg = TRBDF2(autodiff = true)   # fails for fixed dt = 0.01
# alg = Rosenbrock23()

prob = ODEProblem(f_ord, y0, (0.01, 1.0e4))
@time sol = solve(prob, alg, dt = 0.01, reltol = 1e-6,
                #   adaptive = false
                )

if sol.retcode == ReturnCode.Success
    display(plot!(sol.t,  mapreduce(permutedims, vcat, sol.u), 
                  xaxis = :log, color = :black, linewidth = 2, line = :dash)
            )
    @show sol.destats
else 
    @show sol.retcode
end

GC.gc()
println("\ndone")

