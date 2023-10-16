# RKM.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://mjmcnelis.github.io/RKM.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://mjmcnelis.github.io/RKM.jl/dev)

```julia
using RKM
using Plots; plotly()

# differential equation
const C = 0.5
function dy_dt!(f, y; kwargs...)
    f[1] = (y[1] + C) * (1.0 - C - y[1])
    nothing
end

# initial conditions
t0 = -10.0
y0 = exp(t0)/(1.0 + exp(t0)) - C

# final time, initial time step
tf = 10.0
dt0 = 1e-4

# solver options
parameters = Parameters(; method = RungeKutta4(), adaptive = Fixed())

# evolve system
sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, parameters; precision = Float64)

# plot solution
y, t = get_solution(sol)
plot(t, y)
```