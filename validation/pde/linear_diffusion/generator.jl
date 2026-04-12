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
Nx = 20001
x = range(-2000, 2000, Nx)      # grid points
# Nx = 21
# x = range(-10, 10, Nx)      # grid points

dx = x[2] - x[1]            # uniform spacing
p = [a]                     # parameters
abstract_params = dx        # non-sensitivity parameters
# p = [a, dx]               # works for both sensitivity parameters
# abstract_params = nothing

# 0.000259 seconds (32 allocations: 626.172 KiB)
# 0.000068 seconds
# 0.000020 seconds (66 allocations: 177.719 KiB)
# 0.000332 seconds (153 allocations: 357.156 KiB)
# 0.100797 seconds (39.53 k allocations: 115.020 MiB, 2.84% gc time)

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

# plt = plot(sol.t, S, legend = :outertopright, size = (900,500),);
# plt = plot(sol.t[t_idxs], S[t_idxs,:], legend = :outertopright, size = (900,500),);
# plot!(sol.t[t_idxs], SG, label = "", color = :black, line = :dash);
# display(plt)

n_list = [11, 21, 31, 51, 101, 201, 501, 1001, 2001, 5001, 10001, 20001]

# time in microseconds
LU = [5.0, 7.6, 10, 14, 26, 49, 118, 241, 461, 1465, 3031, 6446]

# dim = 30; 1 iteration, 31 applications
# note: much faster if Nx < 30 (is krylovdim <= Nx enforced?)
GX_30 = [11, 34, 63, 66, 76, 99, 171, 290, 533, 1347, 2857, 5562]

# dim = 20; 1 iterations, 21 applications
GX_20 = [11, 31, 31, 34, 39, 51, 91, 160, 312, 751, 1585, 3087]

# dim = 10; 2 iterations, 22 applications (except for Nx = 11 - 31)
GX_10 = [9, 44, 37, 27, 33, 43, 76, 130, 249, 614, 1263, 3633]

# expv defaults
GX_expv = [NaN, NaN, NaN, NaN, NaN, NaN, NaN, 154, 272, 600, 1194, NaN]

plt = plot(legend = :outertopright,
        #    xlims = (1e1, 1e4), ylims = (1e0, 4e3),
        #    xaxis = :log, yaxis = :log,
          );
plot!(n_list, LU, label = "lu!", color = :blue);
plot!(n_list, GX_30, label = "exponentiate (dim = 30)", color = :red);
plot!(n_list, GX_20, label = "exponentiate (dim = 20)", color = :purple);
plot!(n_list, GX_10, label = "exponentiate (dim = 10)", color = :green);
plot!(n_list, GX_expv, label = "expv (dim = 30)", color = :orange);
# display(plt)

GC.gc()
println("\ndone")
