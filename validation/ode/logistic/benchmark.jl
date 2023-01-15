using Revise, RKM, OrdinaryDiffEq
import StaticArrays: SA
import BenchmarkTools: @benchmark, @btime, mean
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/logistic/equations.jl") : nothing

t0 = -10.0
tf = 10.0
dt0_vect = 10.0.^LinRange(-1, -5, 5)

time_ord = Float64[]
memory_ord = Float64[]
time_rkm = Float64[]
memory_rkm = Float64[]

# static = false
static = true

GC.gc()
for dt0 in dt0_vect  
    @show dt0
    # OrdinaryDiffEq
    f = static ? fp_static : fp 
    y0 = static ? SA[exp(t0)/(1.0 + exp(t0)) - 0.5] : [exp(t0)/(1.0 + exp(t0)) - 0.5]
    prob = ODEProblem(f, y0, (t0, tf))
    ord = @benchmark solve($prob, RK4(), dt = $dt0, adaptive = false)
    push!(time_ord, mean(ord).time/1e9)     # convert from ns to s
    push!(memory_ord, ord.memory/1024^2)    # convert from bytes to MiB
    # RKM
    y0 = exp(t0)/(1.0 + exp(t0)) - 0.5
    ps = Parameters(; adaptive = Fixed(), method = RungeKutta4(), 
                      t_range = TimeRange(; t0, tf, dt0))
    rkm = @benchmark evolve_ode($y0, dy_dt!; parameters = $ps, show_progress = false, 
                                static_array = $static)
    push!(time_rkm, mean(rkm).time/1e9)
    push!(memory_rkm, rkm.memory/1024^2)
    GC.gc()
end
# plot options
plot_kwargs = (title = "Logistic equation", titlefontsize = 16, size = (1000, 600),
                 linewidth = 2, legend = :outertopright, legendfontsize = 12,
                 yguidefontsize = 14, ytickfontsize = 12,xguidefontsize = 14, 
                 xtickfontsize = 12, xaxis = :log, yaxis = :log)

# Plot runtimes vs time step
plt = plot(dt0_vect, time_ord; label = "OrdinaryDiffEq", 
           xlabel = "Time step [s]", ylabel = "Runtime [s]",
           xlims = (1e-5, 1e-1), ylims = (1e-5, 1e1), plot_kwargs...)
plot!(dt0_vect, time_rkm; label = "RKM", plot_kwargs...)
display(plt)

# Plot runtimes vs time step
plt = plot(dt0_vect, memory_ord; label = "OrdinaryDiffEq", 
           xlabel = "Time step [s]", ylabel = "Memory [MiB]",
           xlims = (1e-5, 1e-1), ylims = (1e-2, 1e3),
            plot_kwargs...)
plot!(dt0_vect, memory_rkm; label = "RKM", plot_kwargs...)
display(plt)