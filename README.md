# RKM.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://mjmcnelis.github.io/RKM.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://mjmcnelis.github.io/RKM.jl/dev)

## Overview
`RKM.jl` is an ordinary differential equation (ODE) solver written in Julia. It is based from an original version [RKM](https://github.com/mjmcnelis/RKM) written in Python.

This repository is currently in the development/testing phase, with one alpha-version released. The main features of this package include

- Explicit and diagonal-implicit Runge-Kutta methods
- Jacobian evaluation with finite difference or forward auto-differentiation
- Adaptive time stepping (embedded or step doubling) with a PID controller
- Solver statistics (e.g. runtime, solution size, excess memory/allocations)
- Options for float precision, timer, progress bar and static arrays


## Setup
This package is not registered. To install, clone the repository

    cd <your_path_dir>
    git clone https://github.com/mjmcnelis/RKM.jl.git

and develop the package in a Julia REPL (assumes v1.10.2 or higher):
```julia
using Pkg
Pkg.develop(path = raw"<your_path_dir>/RKM.jl")
```

It is also recommended to install these packages in your base environment:
```julia
] add DoubleFloats Plots
```

## Example
This code example shows how to use the ODE solver.
```julia
using RKM
using DoubleFloats: Double64
using Plots; plotly()

# logistic equation
if !(@isdefined dy_dt!)
    function dy_dt!(f, y; p, kwargs...)
        B = p[1]
        f[1] = (y[1] + B) * (1.0 - B - y[1])
        nothing
    end
end

# initial conditions
B = 0.5
t0 = -10.0 |> BigFloat
y0 = exp(t0)/(1.0 + exp(t0)) - B

# final time, initial time step, model parameters
tf = 10.0
dt0 = 1e-4
p = [B]

# solver options
options = SolverOptions(; method = RungeKutta4(), adaptive = Fixed())

# evolve ode, plot solution
precision = Float64
sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options,
                 precision; model_parameters = p)
y, t = get_solution(sol)
plot(t, y)
```
### External inputs
The solver requires a user-defined function `dy_dt!` that numerically evaluates the ODE system at a given state and time
```math
\frac{d\vec{y}}{dt} = \vec{f}(\vec{y}, t; p)
```
We support the following in-place methods to store the result in a vector variable `f`:
```julia
dy_dt!(f, y; kwargs...)
dy_dt!(f, y; t, kwargs...)
dy_dt!(f, y; p, kwargs...)
dy_dt!(f, y; t, p)
```
The first one is for ODEs that depend only on the state variable(s) `y`. The other three are for ODEs that also depend on time `t` and/or model parameters `p`.

In addition, we need to specify the initial state `y0` (either scalar or vector), the initial and final times `t0` and `tf`, the initial time step `dt0`, and model parameters `p` (if any).

### Solver options
Next, we have to set two of the solver options: the ODE method `method` and the adaptive time step algorithm `adaptive`. In this example, we use the classic RK4 method and a fixed time step.

The solver supports a number of explicit and diagonal-implicit Runge-Kutta methods. You can list all of the available methods by calling
```julia
list_explicit_runge_kutta_methods()
list_implicit_runge_kutta_methods()
```
We can set the time step to be either fixed or adaptive. The latter is based on step doubling or embedded techinques.
```julia
Fixed()                                                 # fixed time step
Doubling(; epsilon = 1e-6, alpha = 1e-6, delta = 1e-6)  # step doubling
Embedded(; epsilon = 1e-6, alpha = 1e-6, delta = 1e-6)  # embedded Runge-Kutta
```
All Runge-Kutta methods are compatible with `Fixed` or `Doubling`, whereas `Embedded` is only compatible with embedded methods. The parameters `epsilon`, `alpha` and `delta` control the relative, absolute and incremental error tolerances, respectively.

### ODE evolution
Finally, we call the function `evolve_ode` to evolve the ODE system and store the numerical solution. An in-place version `evolve_ode!` is also available.
```julia
precision = Float64

sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options,
                 precision; model_parameters = p)

sol = Solution(precision)
evolve_ode!(sol, y0, t0, tf, dt0, dy_dt!, options; model_parameters = p)
```
You can adjust the numerical precision of the solver with the argument `precision` (defaulted to `Float64`). For example, we could have used `Double64` or `BigFloat`.

*Note: `model_parameters` can be omitted if `dy_dt!` does not depend on `p`.*

### Post-processing
The field `sol.t` stores the time series ($t_0, ..., t_n$), and `sol.y` stores the solution set ($\vec{y}_0, ..., \vec{y}_n$) in linear column format. If the state vector $\vec{y}$ is one-dimensional, we can plot the solution with
```julia
plot(sol.t, sol.y)
```
If $\vec{y}$ is multi-dimensional, we have to reshape `sol.y` as a (transposed) matrix
```julia
julia> y, t = get_solution(sol);
julia> y
200002×1 transpose(::Matrix{Float64}) with eltype Float64:
 -0.49995460213129755
  ⋮
  0.49995460667064867
```
before plotting the results
```julia
plot(t, y)
```

## Additional features

### Runtime statistics

After the solver finishes, we can print runtime statistics with the function `get_stats`. After recompiling, we get
```julia
julia> @time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options,
                              precision; model_parameters = p);
  0.046932 seconds (152 allocations: 3.063 MiB)

julia> get_stats(sol)
time steps taken     = 200001
time points saved    = 200002
step rejection rate  = 0.0 %
function evaluations = 800004
jacobian evaluations = 0
evolution runtime    = 0.04675 seconds
solution size        = 3.052 MiB
config memory        = 10.344 KiB
excess memory        = 0 bytes
```
Here, we show the number of time steps saved and the number of times `dy_dt!` was evaluated. We also list several performance metrics: the runtime, the solution size, the configuration memory, and the excess memory allocated during the time evolution loop. In this example, almost all of the allocated memory went towards storing the solution.

### Timer and progress display

We can set a time limit and display a progress bar by passing `timer` and `show_progress` to the solver options:
```julia
options_timer = SolverOptions(; method = RungeKutta4(), adaptive = Fixed(),
                                timer = TimeLimit(; wtime_min = 1), # set timer to 1 minute
                                show_progress = true                # display progress
                             )
```
The solver stops if it exceeds the time limit, but it still saves part of the solution.
```julia
julia> dt0_small = 5e-8;             # trigger timer
julia> sol = evolve_ode(y0, t0, tf, dt0_small, dy_dt!, options_timer,
                        precision; model_parameters = p);
Progress:  70%|█████████████████████████████▍            |  ETA: 0:00:25 ( 0.85  s/it)
┌ Warning: Exceeded time limit of 1.0 minutes (stopping evolve_ode!...)
└ @ RKM ~/Desktop/RKM.jl/src/timer.jl:58
```

### Advanced solver options
