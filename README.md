# RKM.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://mjmcnelis.github.io/RKM.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://mjmcnelis.github.io/RKM.jl/dev)

## Overview
`RKM.jl` is an ordinary differential equation (ODE) solver written in Julia. It is currently in the development/testing phase, with several alpha versions released. The main features of this package include:

- Explicit and diagonal-implicit Runge-Kutta methods
- Linear multistep methods (fixed time step only)
- Adaptive time stepping with a PID controller
- Jacobian evaluation with finite differences or forward auto-differentiation
- First-order sensitivity analysis (fixed time step only)
- Arbitrary float precision and dense output
- Options to set a timer and display a progress bar
- High runtime performance and efficient memory usage


## Setup
This package is not registered. To install, clone the repository

    cd <your_path_dir>
    git clone https://github.com/mjmcnelis/RKM.jl.git

and develop the package in a Julia REPL (assumes v1.11.4 or higher):
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
    function dy_dt!(f, y, t; p, kwargs...)
        f[1] = (y[1] + p[1]) * (1.0 - p[1] - y[1])
        return nothing
    end
end

# sensitivity parameters
p = [0.5]

# initial conditions
t0 = -10.0 |> BigFloat
y0 = [exp(t0)/(1.0 + exp(t0)) - p[1]]

# final time, initial time step
tf = 10.0
dt0 = 1e-4

# solver options
options = SolverOptions(; method = RungeKutta4(), adaptive = Fixed(), precision = Float64)

# evolve ode, plot solution
sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p)
t, y = get_solution(sol)
plot(t, y)
```
<details>
<summary>External inputs</summary>

### External inputs
The solver requires a user-defined function `dy_dt!` that numerically evaluates the ODE system at a given state and time
```math
\frac{d\vec{y}}{dt} = \vec{f}(\vec{y}, t; p)
```
We support the following in-place methods to store the result in a vector variable `f`:
```julia
dy_dt!(f, y, t; kwargs...)
dy_dt!(f, y, t; p, kwargs...)
dy_dt!(f, y, t; abstract_params, kwargs...)
dy_dt!(f, y, t; p, abstract_params)
```
The first one is for ODEs that depend on the state variable(s) `y` and time `t`. The other three are for ODEs that also depend on sensitivity and/or abstract parameters.

In addition, we need to specify the initial state vector `y0`, the initial and final times `t0` and `tf`, the initial time step `dt0`, and parameters `(p, abstract_params)` (if any).

*Note: `p` is a `Vector{Float64}` type, while `abstract_params` represents any object.*
</details>

<details>
<summary>Solver options</summary>

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
</details>

<details>
<summary>ODE evolution</summary>

### ODE evolution
Finally, we call the function `evolve_ode` to evolve the ODE system and store the numerical solution. An in-place version `evolve_ode!` is also available.
```julia
sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p)

sol = Solution(options)
evolve_ode!(sol, y0, t0, tf, dt0, dy_dt!, options, p)
```
You can adjust the numerical precision of the solver by changing `precision` in `options` (defaulted to `Float64`). For example, we could have used `Double64` or `BigFloat`.

*Note: `p` can be omitted if `dy_dt!` does not depend on it.*
</details>

<details>
<summary>Post-processing</summary>

### Post-processing
The field `sol.t` stores the time series ($t_0, ..., t_n$), and `sol.y` stores the solution set ($\vec{y}_0, ..., \vec{y}_n$) in linear column format. If the state vector $\vec{y}$ is one-dimensional, we can plot the solution with
```julia
plot(sol.t, sol.y)
```
If $\vec{y}$ is multi-dimensional, we have to reshape `sol.y` as a (transposed) matrix
```julia
julia> t, y = get_solution(sol);
julia> y
200002×1 transpose(::Matrix{Float64}) with eltype Float64:
 -0.49995460213129755
  ⋮
  0.4999546021312957
```
before plotting the results
```julia
plot(t, y)
```
</details>

## Additional features

<details>
<summary>Runtime statistics</summary>

### Runtime statistics
After the solver finishes, we can print runtime statistics with the function `get_stats`. After recompiling, we get
```julia
julia> @time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p);
  0.023003 seconds (149 allocations: 3.060 MiB)

julia> get_stats(sol)
time steps taken     = 200001
time points saved    = 200002
step rejection rate  = 0.0 %
function evaluations = 800005
jacobian evaluations = 0
evolution runtime    = 0.02292 seconds
solution size        = 3.052 MiB
sensitivity size     = 0 bytes
config memory        = 7.641 KiB
excess memory        = 0 bytes
```
Here, we show the number of time steps saved and the number of times `dy_dt!` was evaluated. We also list several performance metrics: the runtime, the solution size, the configuration memory, and the excess memory allocated during the time evolution loop. In this example, almost all of the allocated memory went towards storing the solution.
</details>

<details>
<summary>Timer and progress bar</summary>

### Timer and progress bar
We can set a timer and display a progress bar by passing `timer` and `show_progress` to the solver options:
```julia
options = SolverOptions(; method = RungeKutta4(), adaptive = Fixed(),
                          timer = TimeLimit(; wtime_min = 1), # 1 minute
                          show_progress = true)
dt0_small = 4e-8  # trigger timer
```
The solver stops if it exceeds the time limit, but it still saves part of the solution.
```julia
julia> sol = evolve_ode(y0, t0, tf, dt0_small, dy_dt!, options, p);
  Progress:  82%|███████████████████████████▉      |  ETA: 0:00:13 ( 0.73  s/it)
       runtime: 00:01:00
   total_steps: 413612356
             t: 6.544494206390023
            dt: 4.0e-8
┌ Warning: Exceeded time limit of 1.0 minutes (stopping evolve_ode!...)
└ @ RKM ~/Desktop/RKM.jl/src/timer.jl:108
```
Some variables are displayed in real-time to determine how fast (or slow) the solver is progressing.
</details>

## Advanced solver options

<details>
<summary>Sensitivity analysis</summary>

### Sensitivity analysis

If the ODE function `dy_dt!` depends on $m$ parameters `p`, we can evolve the sensitivity coefficients $\vec{S}_{j} = {\partial\vec{y}/\partial p_j}$ ($j \in m$) with the solver option `sensitivity`.
```julia
options = SolverOptions(; method = RungeKutta4(), adaptive = Fixed(),
                          sensitivity = DecoupledDirect())

sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p)
```
The field `sol.S` stores the sensitivities in linear column format. After reshaping it, we can plot the sensitivity curves $S = [\vec{S}_1(t) ... \vec{S}_m(t)]$
```julia
t, S = get_sensitivity(sol)
plot(t, S)
```
</details>
