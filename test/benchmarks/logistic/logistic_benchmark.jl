using RKM, OrdinaryDiffEq
import StaticArrays: SA
import BenchmarkTools: @benchmark, mean
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/logistic/equations.jl") : nothing
GC.gc()

t0 = -10.0
tf = 10.0

time_ord = Float64[]
memory_ord = Float64[]
time_rkm = Float64[]
memory_rkm = Float64[]

static = false
# static = true
save_solution = true
# save_solution = false

N = 1
y0 = [exp(t0) / (1.0 + exp(t0)) - get_a(i,N) for i in 1:N]
dt0_vect = 10.0.^(-(1:1:5))

# vary time step (fix number of state variables)
for dt0 in dt0_vect
    @show dt0

    # OrdinaryDiffEq
    f = static ? f_ord_static : f_ord
    y = static ? SA[y0...] : y0
    prob = ODEProblem(f, y, (t0, tf))
    saveat = save_solution ? dt0 : ()
    ord = @benchmark solve($prob, RK4(), dt = $dt0, adaptive = false,
                           maxiters = $(10^7), saveat = $saveat,
                        #    save_everystep = $save_solution,
                        #    save_start = $save_solution,
                        #    save_end = $save_solution
                          )
    push!(time_ord, mean(ord).time/1e9)     # convert from ns to s
    push!(memory_ord, ord.memory/1024^2)    # convert from bytes to MiB
    GC.gc()
    # RKM
    options = SolverOptions(; method = RungeKutta4(), adaptive = Fixed(), save_solution)
    rkm = @benchmark evolve_ode($y0, $t0, $tf, $dt0, dy_dt!, $options)
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
           xlims = (1e-5, 1e-1), ylims = (1e-5, 1e1), plot_kwargs...);
plot!(dt0_vect, time_rkm; label = "RKM", plot_kwargs...);
display(plt)

# Plot memory usage vs time step
plt = plot(dt0_vect, memory_ord; label = "OrdinaryDiffEq",
           xlabel = "Time step [s]", ylabel = "Memory usage [MiB]",
           xlims = (1e-5, 1e-1), ylims = (1e-3, 1e3),
            plot_kwargs...);
plot!(dt0_vect, memory_rkm; label = "RKM", plot_kwargs...);
display(plt)

#-----------------------------------------------------------------

time_ord = Float64[]
memory_ord = Float64[]
time_rkm = Float64[]
memory_rkm = Float64[]

dt0 = 1e-1
N_vect = (10).^(0:1:6)

# vary state variables (fix time step)
for N in N_vect
    @show N
    # don't use static arrays
    y0 = [exp(t0) / (1.0 + exp(t0)) - get_a(i,N) for i in 1:N]

    # OrdinaryDiffEq
    prob = ODEProblem(f_ord, y0, (t0, tf))
    saveat = save_solution ? dt0 : ()
    ord = @benchmark solve($prob, RK4(), dt = $dt0, adaptive = false, saveat = $saveat,
                        #    save_everystep = $save_solution,
                        #    save_start = $save_solution,
                        #    save_end = $save_solution
                        )
    push!(time_ord, mean(ord).time/1e9)     # convert from ns to s
    push!(memory_ord, ord.memory/1024^2)    # convert from bytes to MiB
    GC.gc()
    # RKM
    options = SolverOptions(; method = RungeKutta4(), adaptive = Fixed(), save_solution)
    rkm = @benchmark evolve_ode($y0, $t0, $tf, $dt0, dy_dt!, $options)
    push!(time_rkm, mean(rkm).time/1e9)
    push!(memory_rkm, rkm.memory/1024^2)
    GC.gc()
end

# Plot runtimes vs number of state variables
plt = plot(N_vect, time_ord; label = "OrdinaryDiffEq",
           xlabel = "Number of state variables", ylabel = "Runtime [s]",
           xlims = (1e0, 1e6), ylims = (1e-5, 1e2), plot_kwargs...);
plot!(N_vect, time_rkm; label = "RKM", plot_kwargs...);
display(plt)

# Plot memory usage vs number of state variables
plt = plot(N_vect, memory_ord; label = "OrdinaryDiffEq",
           xlabel = "Number of state variables", ylabel = "Memory usage [MiB]",
           xlims = (1e0, 1e6), ylims = (1e-3, 1e4), plot_kwargs...);
plot!(N_vect, memory_rkm; label = "RKM", plot_kwargs...);
display(plt)
