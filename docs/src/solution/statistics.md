
maybe this should be a separate section...

# Solver statistics

After running the ODE solver, you can print out statistics with the function `get_stats`:

```julia
using RKM
using Plots; plotly()

function dy_dt!(f, y, t; p, kwargs...)
    γ = p[1]
    ω = p[2]
    f[1] = y[2]
    f[2] = -γ*y[2] - ω^2*y[1]
    return nothing
end

y0 = [1.0, -1.0]  # [position, velocity]
t0 = 0.0
tf = 10.0
dt0 = 1e-2

γ = 101.0         # damping coefficient
ω = 10.0          # frequency
p = [γ, ω]

options = SolverOptions(; method = RungeKutta4(), adaptive = Fixed(),);

sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p);
@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p);    # recompile
#  0.000208 seconds (162 allocations: 34.250 KiB)
```

```julia
julia> get_stats(sol)
time steps taken     = 1001
time points saved    = 1002
step rejection rate  = 0.0 %
function evaluations = 4005
jacobian evaluations = 0
evolution runtime    = 0.0001559 seconds
solution size        = 23.484 KiB
sensitivity size     = 0 bytes
config memory        = 8.516 KiB
excess memory        = 0 bytes
```

Here we summarize each runtime statistic:

- *Time steps taken*
    - the number of steps the solver advanced with during the time evolution.
- *Time points saved*
    - the number of solution points saved if you chose to output them (i.e. `save_solution = true`)
- *Step rejection rate*
    - the percentage of attempted steps that were rejected if you used an adaptive time step method.
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
    - The memory allocated to configure the solver before running the time evolution loop (mostly for intermediate caches).
- *Excess memory*
    - The amount of memory allocated during the time evolution loop (discounting solution storage).
