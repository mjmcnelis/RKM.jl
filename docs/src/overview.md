# Overview

```julia
using RKM
import DoubleFloat
using Plots; plotly()

# differential equation
const C = 0.5
function dy_dt!(f, t, y)
    f[1] = (y[1] + C) * (1.0 - C - y[1])
    nothing
end

# initial conditions
t0 = -10.0
y0 = exp(t0)/(1.0 + exp(t0)) - C

# time range, solver options
t_range = TimeRange(; t0, tf = 10.0, dt0 = 1e-4)
options = SolverOptions(; t_range, method = RungeKutta4(), adaptive = Fixed())

# evolve system
sol = evolve_ode(y0, dy_dt!, options, precision = Float64)

# plot solution
t, y = get_solution(sol)
plot(t, y)
```