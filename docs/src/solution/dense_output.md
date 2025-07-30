
# Dense output

Interpolation can only be performed after solving the ODE, since it uses raw solution data as input.
In addition to the state variables, the cubic Hermite and continuous formula routines use time derivative and intermediate stage data, respectively. You will need to change the solver option `interpolator` to solve and store the solution variables required by each routine.

## No interpolation

If you do not plan to interpolate the solution, use the default value `interpolator = NoInterpolator()`. The solver will only output the time series and state variables of the raw solution.

## Cubic Hermite

All ODE methods can use cubic Hermite interpolation to generate dense output. Here, we continue with the overdamped oscillator example and set the solver option `interpolator = CubicHermite()`. The solver will output the time derivatives, which is stored in `sol.f`.

```julia
options = SolverOptions(; method = TrapezoidRuleBDF21(),
                          adaptive = Embedded(; alpha = 1e-3),
                          interpolator = CubicHermite(),);

sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p);
```

First, we plot the non-uniform solution data from the adaptive TRBDF2 method:

```julia
t, y = get_solution(sol)
scatter(t, y; xlabel = "t", ylabel = "y", label = ["x" "v"], color = [:red :blue], ms = 3)
```

```@raw html
<img src="scatter_plot.png" width="600">
```

Next, we call `interpolate_solution` to interpolate the data uniformly (the keyword argument `dt_dense` controls the temporal spacing).

```julia
t_dense, y_dense = interpolate_solution(options, sol; dt_dense = 1e-4);
plot!(t_dense, y_dense; label = ["x (dense)" "v (dense)"], color = [:red :blue])
```

The interpolated solution is $C^1$ continuous and third-order accurate (i.e. it has a local error of $\mathcal{O}(\Delta t_n^4)$ between the raw time intervals $[t_n, t_{n+1}]$).

```@raw html
<img src="hermite_plot.png" width="600">
```

## Continuous formula

```julia
options = SolverOptions(; method = TrapezoidRuleBDF21(),
                          adaptive = Embedded(; alpha = 1e-3),
                          interpolator = ContinuousFormula(),);

sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p);

t, y = get_solution(sol)
scatter(t, y; xlabel = "t", ylabel = "y", label = ["x" "v"], color = [:red :blue], ms = 3)
```

```julia
t_dense, y_dense = interpolate_solution(options, sol; dt_dense = 1e-4);
plot!(t_dense, y_dense; label = ["x (dense)" "v (dense)"], color = [:red :blue])
```