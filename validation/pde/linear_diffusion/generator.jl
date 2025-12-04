using Revise, RKM, LinearSolve, ForwardDiff
using AppleAccelerate
using DoubleFloats
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/pde/linear_diffusion/equations_dirichlet.jl") : nothing

precision = Float64
# precision = Double64

a = 0.25                    # diffusion constant

# 0.004267 seconds (8.81 k allocations: 4.499 MiB)

# note: if try to go to larger matrices, the CFL can be bad
# Nx = 1001
# x = range(-100, 100, Nx)      # grid points
Nx = 21
x = range(-10, 10, Nx)      # grid points

dx = x[2] - x[1]            # uniform spacing
p = [a]                     # parameters
abstract_params = dx        # non-sensitivity parameters
# p = [a, dx]               # works for both sensitivity parameters
# abstract_params = nothing

t0 = 1.0                    # initial conditions
y0 = gauss.(x, t0; p, t0)

tf = 100.0
dt0 = 0.01

CFL = 2.0*a*dt0/dx^2        # CFL number
@show CFL
Nt = 300                    # temporal stride (for plot)

# generate sparsity pattern via nansafe
sparsity = nansafe_state_jacobian(y0, t0, dy_dt!, p; chunk_size = 1, abstract_params)

if precision == Float64
    method = BackwardEuler1()
else
    method = Ketcheson4()
end

options = SolverOptions(;
              method = method,
              adaptive = Fixed(),
              state_jacobian = ForwardColorJacobian(; sparsity),
              root_finder = Newton(; linear_method = KLUFactorization(),),
              sensitivity = DecoupledDirect(; param_jacobian = ForwardJacobian(),),
              time_subroutine = true,
              precision,)

@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p; abstract_params)
_, y = get_solution(sol)
_, S = get_sensitivity(sol)

# get_subroutine_times(sol)

A = sparsity    # TODO: use state_jacobian w/ (t0, y0)
# A = Matrix(sparsity)
@time SG, t_idxs = post_generator(sol, options, dy_dt!, A, p; abstract_params)

plt = plot(sol.t, S, legend = :outertopright, size = (900,500),);
# plt = plot(sol.t[t_idxs], S[t_idxs,:], legend = :outertopright, size = (900,500),);
plot!(sol.t[t_idxs], SG, label = "", color = :black, line = :dash);
# display(plt)

GC.gc()
println("\ndone")
