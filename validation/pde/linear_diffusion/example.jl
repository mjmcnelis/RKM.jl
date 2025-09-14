using Revise, RKM, LinearSolve, ForwardDiff
using AppleAccelerate
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/pde/linear_diffusion/equations.jl") : nothing

show_plot = false            # plot solution
benchmark_diffeq = false    # compare to OrdinaryDiffEq

a = 0.25                    # diffusion constant
Nx = 201
x = range(-10, 10, Nx)      # grid points
dx = x[2] - x[1]            # uniform spacing
p = [a, dx]                 # parameters

t0 = 1.0                    # initial conditions
y0 = gauss.(x, t0; p, t0)

tf = 10.0
dt0 = 0.01

CFL = 2.0*a*dt0/dx^2        # CFL number
@show CFL
Nt = 300                    # temporal stride (for plot)

# generate sparsity pattern via nansafe
sparsity = nansafe_state_jacobian(y0, t0, dy_dt!, p; chunk_size = 1)
display(sparsity)

# test_nansafe(; x = ForwardDiff.Dual(NaN, 1.0))
# test_nansafe(; x = ForwardDiff.Dual(NaN, 0.0))
# q()

options = SolverOptions(
              method = BackwardEuler1(),
              adaptive = Fixed(),
              # TODO: finitediff w/ sparsity doesn't seem to be working
            #   state_jacobian = FiniteJacobian(; sparsity),
              state_jacobian = ForwardJacobian(; sparsity),
              root_finder = Newton(; linear_method = KLUFactorization(),),
              time_subroutine = true
          )

@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p)
get_stats(sol)
get_subroutine_times(sol)

GC.gc()

# OrdinaryDiffEq
if benchmark_diffeq
    println("")
    using OrdinaryDiffEq

    # prob = ODEProblem(diffeq!, y0, (t0, tf), p)
    func = ODEFunction(diffeq!; jac_prototype = sparsity)
    prob = ODEProblem(func, y0, (t0, tf), p)

    alg = ImplicitEuler(autodiff = true, linsolve = KLUFactorization())

    @time sol_diffeq = solve(prob, alg, dt = dt0, adaptive = false);
    # display(sol_diffeq.destats)
    @show sol_diffeq.destats.nf sol_diffeq.destats.njacs
end

if show_plot
    # initial time slice
    plt = plot(x, y0, label = "t = $t0", color = "indianred", linewidth = 2,
            size = (900, 600), ylims = (-0.5, 1.3),
            ylabel = "u", yguidefontsize = 14, ytickfontsize = 12,
            xlabel = "t", xguidefontsize = 14, xtickfontsize = 12,
            legend = :outertopright, legendfontsize = 12)

    t, y = get_solution(sol)
    for n = 1:3
        t_slice = t[1 + Nt*n]
        plot!(x, y[1 + Nt*n, :], color = "indianred", linewidth = 2, label = "t = $t_slice")
    end

    plot!(x, y0, label = "t = $t0 (exact)", color = "black", linewidth = 1, line = :dash)
    for n = 1:3
        t_slice = t[1 + Nt*n]
        y_exact = gauss.(x, t_slice; p, t0)
        plot!(x, y_exact, color = "black", linewidth = 1,
              line = :dash, label = "t = $t_slice (exact)")
    end
    display(plt)
end

GC.gc()
println("\ndone")

