using Revise, RKM, OrdinaryDiffEq
import StaticArrays: SA
import BenchmarkTools: @benchmark, @btime, mean
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/gaussian/equations.jl") : nothing
GC.gc()

t0 = -5.0
y0 = [exp(-t0^2/2.0) + 1.0]
dt0 = 1e-4
tf = 5.0

safety = 0.9
low = 0.1
high = 5.0

time_ord = Float64[]
memory_ord = Float64[]
time_rkm = Float64[]
memory_rkm = Float64[]

static = false
# static = true

epsilon_vect = 10.0.^(-(10:5:15))

# vary time step (fix number of state variables)
for epsilon in epsilon_vect
    @show epsilon

    # OrdinaryDiffEq
    f = static ? f_ord_static : f_ord
    y = static ? SA[y0...] : y0
    prob = ODEProblem(f, y, (t0, tf))

    ord = @benchmark solve($prob, DP5(), dt = $dt0, reltol = $epsilon, abstol = 0.0,
                           qmin = 1/$high, qmax = 1/$low, gamma = $safety,
                           controller = IController())

    push!(time_ord, mean(ord).time/1e9)     # convert from ns to s
    push!(memory_ord, ord.memory/1024^2)    # convert from bytes to MiB
    GC.gc()
    # RKM
    ps = Parameters(; adaptive = Embedded(; epsilon, low, high, safety),
                      method = DormandPrince54(), t_range = TimeRange(; t0, tf, dt0))

    rkm = @benchmark evolve_ode($y0, dy_dt!; parameters = $ps, show_progress = false,
                                static_array = $static)
    push!(time_rkm, mean(rkm).time/1e9)
    push!(memory_rkm, rkm.memory/1024^2)
    GC.gc()
end

# plot options
plot_kwargs = (title = "Gaussian equation", titlefontsize = 16, size = (1000, 600),
                 linewidth = 2, legend = :outertopright, legendfontsize = 12,
                 yguidefontsize = 14, ytickfontsize = 12,xguidefontsize = 14,
                 xtickfontsize = 12, xaxis = :log, yaxis = :log)

# Plot runtimes vs time step
plt = plot(epsilon_vect, time_ord; label = "OrdinaryDiffEq",
           xlabel = "Relative tolerance", ylabel = "Runtime [s]",
           xlims = (1e-20, 1e-10), ylims = (1e-5, 1e1), plot_kwargs...);
plot!(epsilon_vect, time_rkm; label = "RKM", plot_kwargs...);
display(plt)

# Plot memory usage vs time step
plt = plot(epsilon_vect, memory_ord; label = "OrdinaryDiffEq",
           xlabel = "Relative tolerance", ylabel = "Memory usage [MiB]",
           xlims = (1e-20, 1e-10), ylims = (1e-2, 1e3),
            plot_kwargs...);
plot!(epsilon_vect, memory_rkm; label = "RKM", plot_kwargs...);
display(plt)
