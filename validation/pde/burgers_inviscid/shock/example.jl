using Revise, RKM, LinearSolve
using FiniteDiff: finite_difference_jacobian!
using SparseArrays: sparse
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/pde/burgers_inviscid/shock/equations.jl") : nothing

# dy_dt! = dy_dt_central!
# dy_dt! = dy_dt_lax_friedrichs!
# dy_dt! = dy_dt_murman_roe!
dy_dt! = dy_dt_rusanov!

show_plot = true            # plot solution

Nx = 101
x = range(-5, 5, Nx)        # grid points
dx = x[2] - x[1]            # uniform spacing
dt = 0.05                   # time step
p = [dx, dt]                # parameters

t0 = 0.0                    # initial conditions
y0 = shock.(x)

tf = 6.0
dt0 = dt

CFL = dt/(2dx)              # CFL number, constant speed = (Δy/2)
@show CFL

Nt = 40                     # temporal stride (for plot)

options = SolverOptions(; adaptive = Fixed(),
                          method = Euler1(),
                        #   method = BackwardEuler1(),
                          interpolator = CubicHermite()#; dt_save = 0.05)
                       )

@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p)
t, y = get_solution(sol)

get_stats(sol)

plt = plot(x, y0, label = "t = 0", color = "indianred", linewidth = 2,
           size = (900, 600), ylims = (-0.25, 1.25),
           ylabel = "u", yguidefontsize = 14, ytickfontsize = 12,
           xlabel = "x", xguidefontsize = 14, xtickfontsize = 12,
           legend = :outertopright, legendfontsize = 12, dpi = 200)

for i = 1:3
    t_slice = t[1 + Nt*i]
    plot!(x, y[1 + Nt*i, :], color = "indianred", linewidth = 2, label = "t = $t_slice")
end

plot!(x, y0, label = "t = 0 (exact)", color = "black", linewidth = 1, line = :dash)
for i = 1:3
    t_slice = t[1 + Nt*i]
    y_exact = shock.(x .- t_slice/2.0)
    plot!(x, y_exact, color = "black", linewidth = 1,
          line = :dash, label = "t = $t_slice (exact)")
end
display(plt)

println("\ndone")
