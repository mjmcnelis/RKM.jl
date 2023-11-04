# RKM.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://mjmcnelis.github.io/RKM.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://mjmcnelis.github.io/RKM.jl/dev)

## Overview
RKM.jl is an ordinary differential equation (ODE) solver written in the Julia language. It is based off an earlier version [RKM](https://github.com/mjmcnelis/RKM) written in Python.

This repository is currently under development, with one alpha-version released. The package contains the following main features:

- Explicit and diagonal-implicit Runge-Kutta methods
- Evaluate Jacobian with finite difference or forward auto-differentiation
- Adaptive time stepping (embedded or step doubling) with a PID controller
- Options for timer, progress bar, float precision and static arrays
- Solver statistics: runtime, solution size and excess memory/allocations

## Setup
This package is currently not registered. To install, clone the repository

    cd <your_path_dir>
    git clone https://github.com/mjmcnelis/RKM.jl.git

and do the following in a Julia REPL (assumes v1.9.0 or higher):
```julia
using Pkg
Pkg.develop(path = raw"<your_path_dir>/RKM.jl")
```

It is also recommended to install these external packages:
```julia
] add Plots DoubleFloats
```

## Example
The following is a basic code example of the ODE solver's usage:
```julia
using RKM
using Plots; plotly()

# differential equation
function dy_dt!(f, y; p, kwargs...)
    C = p[1]
    f[1] = (y[1] + C) * (1.0 - C - y[1])
    nothing
end

# initial conditions, model parameters
C = 0.5
t0 = -10.0
y0 = exp(t0)/(1.0 + exp(t0)) - C
p = [C]

# final time, initial time step
tf = 10.0
dt0 = 1e-4

# solver options
parameters = Parameters(; method = RungeKutta4(), adaptive = Fixed())

# evolve system
sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, parameters;
                 model_parameters = p, precision = Float64)

# plot solution
y, t = get_solution(sol)
plot(t, y)
```
### External inputs
The package requires you to define a function `dy_dt!` that numerically evaluates the ODE system at a given state/time
```math
\frac{dy_i}{dt} = f_i(\vec{y}; t, p)
```
We support the following in-place methods to store the result in a vector variable `f`:
```julia
dy_dt!(f, y; kwargs...)
dy_dt!(f, y; t, kwargs...)
dy_dt!(f, y; p, kwargs...)
dy_dt!(f, y; t, p)
```
The first method is for ODEs that depend only on the state variable(s) `y`. You can pass in keyword arguments if your ODE system also depends on time `t` and/or model parameters `p`.

In addition, you will need to specify the initial conditions `y0` (either scalar or vector), the initial and final times `t0` and `tf`, the initial time step `dt0` and model parameters `p` (if any).