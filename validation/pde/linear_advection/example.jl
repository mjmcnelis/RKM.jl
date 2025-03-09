using Revise, RKM, LinearSolve
using FiniteDiff: finite_difference_jacobian!
using SparseArrays: sparse
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/pde/linear_advection/equations.jl") : nothing

# dy_dt! = dy_dt_central!
# dy_dt! = dy_dt_lax_friedichs!
# dy_dt! = dy_dt_upwind!
dy_dt! = dy_dt_rusanov!

show_plot = true           # plot solution

a = 1.0                     # wave speed
Nx = 121
x = range(0, 12, Nx)        # grid points
dx = x[2] - x[1]            # uniform spacing
dt = 0.05                   # time step
p = [a, dx, dt]             # parameters

t0 = 0.0                    # initial conditions
y0 = gauss.(x)

tf = 6.0
dt0 = dt

CFL = a*dt/dx               # CFL number
@show CFL

Nt = 40                     # temporal stride (for plot)

options = SolverOptions(; # note: LxF is only compatible w/ constant dt
                          adaptive = Fixed(),
                          method = Euler1(),
                        #   method = BackwardEuler1(),
                          interpolator = CubicHermite()#; dt_save = 0.05)
                       )

@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p)
t, y = get_solution(sol)

get_stats(sol)

if show_plot
    plt = plot(x, y0, label = "t = 0", color = "indianred", linewidth = 2,
            size = (900, 600), ylims = (-0.5, 1.3),
            ylabel = "u", yguidefontsize = 14, ytickfontsize = 12,
            xlabel = "x", xguidefontsize = 14, xtickfontsize = 12,
            legend = :outertopright, legendfontsize = 12, dpi = 200)
    for n = 1:3
        t_slice = t[1 + Nt*n]
        plot!(x, y[1 + Nt*n, :], color = "indianred", linewidth = 2, label = "t = $t_slice")
    end

    plot!(x, y0, label = "t = 0 (exact)", color = "black", linewidth = 1, line = :dash)
    for n = 1:3
        t_slice = t[1 + Nt*n]
        y_exact = gauss.(x .- a*t_slice)
        plot!(x, y_exact, color = "black", linewidth = 1,
            line = :dash, label = "t = $t_slice (exact)")
    end
    display(plt)
end

println("\ndone")
