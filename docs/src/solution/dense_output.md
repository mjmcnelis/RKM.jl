
# Dense output

The interpolation routines depend on raw solution data as input. You will need to adjust the default solver option `interpolator = NoInterpolator()` to solve and store the solution variables required by each routine.

## Cubic Hermite

All ODE methods can use cubic Hermite interpolation to generate dense output. Here, we continue with the overdamped oscillator example and set the solver option `interpolator = CubicHermite()`.

```julia
options = SolverOptions(; method = TrapezoidRuleBDF2(), adaptive = Doubling(),
                          interpolator = CubicHermite(),);

sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p);
```

We also use the TRBDF2 method with step-doubling to produce non-uniform solution data that is third-order accurate:

```julia
t, y = get_solution(sol)
scatter(t, y; xlabel = "t", ylabel = "y", label = ["x" "v"], color = [:red :blue], ms = 3)
```

```@raw html
<img src="scatter_plot.png" width="600">
```

The solution should contain state variable and time derivative data, which are needed for cubic Hermite interpolation. Now we can call `interpolate_solution` to interpolate the data uniformly (the keyword argument `dt_dense` changes the temporal spacing).

```julia
t_dense, y_dense = interpolate_solution(options, sol; dt_dense = 1e-4);
plot!(t_dense, y_dense; label = ["x (dense)" "v (dense)"], color = [:red :blue])
```

```@raw html
<img src="hermite_plot.png" width="600">
```

## Continuous formula
