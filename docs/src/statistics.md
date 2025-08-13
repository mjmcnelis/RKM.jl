
## Solver statistics

After evolving an ODE system, you can print out runtime statistics to analyze the solver's performance. In this example, we solve the 1D diffusion equation with central differences:

```math
    \frac{dT_i}{dt} = \frac{\alpha(T_{i+1} - 2T_i + T_{i-1})}{\Delta x^2}
```
where $T_i$ is the temperature at grid point $i$, $\alpha$ is the thermal diffusivity and $\Delta x$ is the grid spacing.

```julia
using RKM

function dy_dt!(f, y, t; p, kwargs...)
    a = p[1]
    dx = p[2]
    N = length(y)
    for i in 1:N
        im1 = max(i-1, 1)  # Neumann boundary conditions
        ip1 = min(i+1, N)
        f[i] = a * (y[ip1] - 2.0*y[i] + y[im1]) / dx^2
    end
    return nothing
end

# parameters
a = 0.25                    # thermal diffusivity
dx = 0.1                    # grid spacing
p = [a, dx]

# initial conditions
t0 = 1.0
x = -10.0:dx:10.0           # 201 grid points
y0 = exp.(-x.^2.0./(4.0*a*t0));

tf = 10.0
dt0 = 0.01

options = SolverOptions(; method = BackwardEuler1(), adaptive = Fixed(),);

sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p);
GC.gc() # garbage collection
```

After recompiling, you can call the function `get_stats` to print the statistics.

```julia
julia> @time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p);
  0.277305 seconds (4.68 k allocations: 284.812 MiB, 1.68% gc time)

julia> get_stats(sol)
time steps taken     = 901
time points saved    = 902
step rejection rate  = 0.0 %
function evaluations = 183805
jacobian evaluations = 901
evolution runtime    = 0.2772 seconds
solution size        = 1.390 MiB
sensitivity size     = 0 bytes
config memory        = 344.336 KiB
excess memory        = 283.061 MiB
```

Here we summarize each runtime statistic:

- *Time steps taken*
    - the number of steps the solver advanced with during the time evolution.
- *Time points saved*
    - the number of solution points saved if you chose to output them (i.e. `save_solution = true`)
- *Step rejection rate*
    - the percentage of attempted steps that were rejected if you used an adaptive time step.
- *Function evaluations*
    - the number of times the ODE function $\frac{d\vec{y}}{dt} = \vec{f}(\vec{y}, t; \vec{p})$ was called. This includes any function evaluations needed the compute the state and parameter Jacobians $\frac{\partial f_i}{\partial y_j}$ and $\frac{\partial f_i}{\partial p_j}$.
- *Jacobian evaluations*
    - the number of times the state and parameter Jacobians were evaluated if you used an implicit solver with Newton's method and/or a sensitivity method.
- *Evolution runtime*
    - the portion of runtime the solver took running the time evolution loop.
- *Solution size*
    - the amount of memory used to store the solution `sol.t` and `sol.y`. If you used an interpolator, the time derivatives `sol.f` or intermediate stages `sol.dy` are also counted.
- *Sensitivity size*
    - the amount of memory used to store the sensitivity coefficients `sol.S` if you used a sensitivity method.
- *Configuration memory*
    - the memory allocated to configure the solver before running the time evolution loop (mostly for intermediate caches).
- *Excess memory*
    - the amount of memory allocated during the time evolution loop (discounting solution storage).

## Excess allocations

We designed the solver to minimize excess allocations during the time evolution. However, there are a few scenarios where excess memory will register:

- the ODE function `dy_dt!` provided by the user can allocate; this can provide insight on how to optimize it.

- in implicit solvers, the default Newton method allocates for each LU factorization of the state Jacobian.

- if you use an adaptive time step, the solution data is resized during each step.

In this example, all of the excess memory results from the second scenario.

## Subroutine runtimes

You can get runtime estimates for several core routines if you set the solver option `benchmark_subroutines = true`.

```julia
options = SolverOptions(; method = BackwardEuler1(), adaptive = Fixed(),
                          benchmark_subroutines = true,);

sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p);
```

After recompiling, we get

```julia
julia> @time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p);

  Subroutine times (seconds)
---------------------------------
function evaluations | 0.06058
jacobian evaluations | 0.04526
linear solve         | 0.2056
save solution        | 0.0001892

  0.299631 seconds (5.20 k allocations: 313.091 MiB, 3.01% gc time)
```

We can see that most of the computational time is spent solving linear systems for the Newton iterations (linear solve times are $\mathcal{O}(n_y^3)$ since the Jacobian matrix is dense by default). This suggests that we should try using a sparse linear solver to reduce the runtime. Following section **(X)**, we generate a sparsity pattern with the function `nansafe_state_jacobian`.

*Note: this assumes that you already set* `NANSAFE_MODE_ENABLED = true` *in* `ForwardDiff`.

```julia
julia> sparsity = nansafe_state_jacobian(y0, t0, dy_dt!, p; chunk_size = 1)
201×201 SparseArrays.SparseMatrixCSC{Float64, Int64} with 601 stored entries:
⎡⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎤
⎢⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⠀⠀⎥
⎢⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⡀⠀⎥
⎣⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠈⠻⣦⎦
```

```julia
using LinearSolve

options = SolverOptions(; method = BackwardEuler1(), adaptive = Fixed(),
                          state_jacobian = FiniteJacobian(; sparsity),
                          root_finder = Newton(; linear_method = KLUFactorization(),),
                          benchmark_subroutines = true,);

@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p);
```

Recompiling the solver yields