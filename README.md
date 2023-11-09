# RKM.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://mjmcnelis.github.io/RKM.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://mjmcnelis.github.io/RKM.jl/dev)

## Overview
`RKM.jl` is an ordinary differential equation (ODE) solver written in Julia. It is based from an original version [RKM](https://github.com/mjmcnelis/RKM) written in Python.

This repository is currently in the development/testing phase, with one alpha-version released. The main features of this package include

- Explicit and diagonal-implicit Runge-Kutta methods
- Evaluate Jacobian with finite difference or forward auto-differentiation
- Adaptive time stepping (embedded or step doubling) with a PID controller
- Options for float precision, static arrays, progress bar and timer
- Solver statistics (e.g. runtime, solution size, excess memory/allocations)

## Setup
This package is not registered. To install, clone the repository

    cd <your_path_dir>
    git clone https://github.com/mjmcnelis/RKM.jl.git

and develop the package in a Julia REPL (assumes v1.9.0 or higher):
```julia
using Pkg
Pkg.develop(path = raw"<your_path_dir>/RKM.jl")
```

It is also recommended to install these packages in your base environment:
```julia
] add Plots DoubleFloats
```

## Example
This code example shows how to use the ODE solver.
```julia
using RKM
using Plots; plotly()
using DoubleFloats: Double64

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
parameters = Parameters(; method = RungeKutta4(), adaptive = Fixed())

# evolve system, print stats
@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, parameters;
                       model_parameters = p, precision = Float64)
get_stats(sol)

# plot solution
@time y, t = get_solution(sol)
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
The first one is for ODEs that depend only on the state variable(s) `y`. If the ODE system also depends on time `t` and/or model parameters `p`, they can be passed as keyword arguments.

In addition, we need to specify the initial state `y0` (either scalar or vector), the initial and final times `t0` and `tf`, the initial time step `dt0`, and model parameters `p` (if any).

### Solver options
Next, we have to set two of the solver options: the ODE method `method` and the adaptive time step `adaptive`. In this example, we use the classic RK4 method and a fixed time step.

The solver supports a number of explicit and diagonal-implicit Runge-Kutta methods. You can list all of the available methods by calling
```julia
list_explicit_runge_kutta_methods()
list_implicit_runge_kutta_methods()
```

We can set the time step to be either fixed or adaptive. The latter is based on step doubling or embedded techinques. All Runge-Kutta methods are compatible with `adaptive = Fixed()` or `Doubling()`, whereas `adaptive = Embedded()` is only compatible with embedded methods.

### ODE evolution
Finally, we call the function `evolve_ode` to evolve the ODE system and store the numerical solution ```y(t)```. An in-place version `evolve_ode!` is also available.
```julia
sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, parameters;
                 model_parameters = p, precision = Float64)

sol = Solution(; precision = Float64)
evolve_ode!(sol, y0, t0, tf, dt0, dy_dt!, parameters; model_parameters = p)
```
We can adjust the numerical precision of the solver with the keyword argument `precision` (e.g. `Float64`, `Double64`, `BigFloat`).

*Note: `model_parameters` can be omitted if `dy_dt!` does not depend on `p`.*

### Runtime statistics
```
# after recompiling
julia> @time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, parameters;
                              model_parameters = p, precision = Float64);
0.017567 seconds (306 allocations: 3.076 MiB)
julia> get_stats(sol)
time steps           = 200002
step rejection rate  = 0.0 %
function evaluations = 800004
jacobian evaluations = 0
solver runtime       = 0.01725 seconds
solution storage     = 3.052 MiB
excess memory        = 0 bytes
excess allocations   = 0
```
