# RKM.jl

```julia 
using RKM 
using Plots; plotly() 

# differential equation
const a = 0.5
function dy_dt!(f, t, y)
    f[1] = (y[1] + a) * (1.0 - a - y[1])
    nothing
end

# initial conditions
t0 = -10.0
y0 = exp(t0) / (1.0 + exp(t0)) - a

# parameters and solver options
adaptive = Fixed()
method = RungeKutta4()
t_span = TimeSpan(; t0, tf = 10.0, dt0 = 1e-4) # time interval
parameters = Parameters(; adaptive, method, t_span)

# evolve ODE
sol = evolve_ode(y0, dy_dt!; parameters)

# plot solution
y, t = get_solution(sol)
plot(t, y)
```