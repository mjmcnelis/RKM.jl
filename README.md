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
This code example shows how to use the ODE solver:
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
model_parameters = [C]

# final time, initial time step
tf = 10.0
dt0 = 1e-4

# solver options
parameters = Parameters(; method = RungeKutta4(), adaptive = Fixed())

# evolve system
sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, parameters;
                 model_parameters, precision = Float64)

# plot solution
y, t = get_solution(sol)
plot(t, y)
```