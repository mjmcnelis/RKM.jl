using Revise, RKM, LinearSolve
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/pde/linear_diffusion/equations.jl") : nothing

show_plot = true            # plot solution

a = 0.25                    # diffusion constant
x = LinRange(-10, 10, 201)  # grid points
dx = x[2] - x[1]            # uniform spacing
p = [a, dx]                 # parameters

# initial conditions
t0 = 1.0
y0 = gauss.(x, t0; p, t0)

tf = 10.0
dt0 = 0.01

CFL = 2.0*a*dt0/dx^2        # CFL number
@show CFL
Nt = 300                    # temporal stride

options = SolverOptions(
              method = BackwardEuler1(),
              adaptive = Fixed(),
              stage_finder = ImplicitStageFinder(
                                 linear_method = LUFactorization(),
                                #  jacobian_method = ForwardJacobian(),
                             ),
          )

@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p)
get_stats(sol)

if show_plot
    # initial time slice
    plt = plot(x, y0, label = "t = $t0", color = "indianred", linewidth = 2,
            size = (900, 600), ylims = (-0.5, 1.3),
            ylabel = "u", yguidefontsize = 14, ytickfontsize = 12,
            xlabel = "t", xguidefontsize = 14, xtickfontsize = 12,
            legend = :outertopright, legendfontsize = 12)

    t, y = get_solution(sol)
    for i = 1:3
        global t = t0 + Nt*i*dt0
        plot!(x, y[1 + Nt*i, :], color = "indianred", linewidth = 2, label = "t = $t")
    end

    plot!(x, y0, label = "t = $t0 (exact)", color = "black", linewidth = 1, line = :dash)
    for i = 1:3
        global t = t0 + Nt*i*dt0
        y_exact = gauss.(x, t; p, t0)
        plot!(x, y_exact, color = "black", linewidth = 1, line = :dash, label = "t = $t (exact)")
    end
    display(plt)
end

GC.gc()
println("\ndone")

