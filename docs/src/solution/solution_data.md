
# Solution data

After solving the ODE, you can plot the time series data stored in the `Solution` struct `sol`.
The following sections are based on the overdamped harmonic oscillator:

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

ω = 10.0          # frequency
γ = 101.0         # damping coefficient
p = [γ, ω]

options = SolverOptions(; method = RungeKutta4(), adaptive = Fixed(),);

sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p);
```

## Time and state variables

The time series $t = (t_0, ..., t_n)$ is stored in `sol.t`. The state variables `sol.y` are stored in linear column format as $\vec{y}(t) = (\vec{y}_0, ..., \vec{y}_n)$.

```julia
julia> sol.t
1002-element Vector{Float64}:
  0.0
  ⋮
 10.0
julia> sol.y
2004-element Vector{Float64}:
  1.0
 -1.0
  ⋮
  4.5399929800627256e-5
 -4.5399929800627256e-5
```

To plot the state variables versus time, we reshape `sol.y` into a transposed matrix by calling the function `get_solution`:

```julia
julia> t, y = get_solution(sol);
julia> y
1002×2 transpose(::Matrix{Float64}) with eltype Float64:
 1.0         -1.0
 ⋮
 4.53999e-5  -4.53999e-5
```

Each column in `y` corresponds to the time series solution of each state variable. Then you can plot the solution with
```julia
plot(t, y; xlabel = "t", ylabel = "y", label = ["x" "v"])
```

```@raw html
<img src="y_plot.png" width="600">
```

## Time derivatives

You can save the first-order time derivatives $\vec{f} = d\vec{y}/dt$ by setting the solver option  `save_time_derivative = true`.
```julia
options = SolverOptions(; method = RungeKutta4(), adaptive = Fixed(),
                          save_time_derivative = true,);

sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p);
```
*Note: using the solver option* `interpolator = CubicHermite()` *will also output the time derivatives.*

The time derivative data is stored in `sol.f` as a linear column (same as `sol.y`). You can plot them by doing

```julia
t, f = get_time_derivative(sol);
plot(t, f; xlabel = "t", ylabel = "f", label = ["dx/dt" "dv/dt"])
```

```@raw html
<img src="f_plot.png" width="600">
```

## Sensitivity coefficients

If you used the solver option `sensitivity = DecoupledDirect()`, you can plot the first-order sensitivity coefficients $\vec{S}_{j} = {\partial\vec{y}/\partial p_j}$.

```julia
options = SolverOptions(; method = RungeKutta4(), adaptive = Fixed(),
                          sensitivity = DecoupledDirect(),);

sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p);
```

The sensitivity coefficients are stored in `sol.S` as a linear column. After reshaping it, you get

```julia
julia> t, S = get_sensitivity(sol);
julia> S
1002×4 transpose(::Matrix{Float64}) with eltype Float64:
 0.0          0.0          0.0         0.0
 ⋮
 4.58122e-6  -4.12263e-6  -9.16244e-5  8.24527e-5
```

The first two columns are the state variables' sensitivity to the first parameter $γ$ (the last two columns are the sensitivity w.r.t. the second parameter $ω$).

```julia
plot(t, S; xlabel = "t", ylabel = "S", label = ["dx/dγ" "dv/dγ" "dx/dω" "dv/dω"])
```

```@raw html
<img src="S_plot.png" width="600">
```

## API Reference

```@autodocs
Modules = [RKM]
Pages   = ["src/solution/solution.jl",
           "src/solution/sizehint.jl"]
```
