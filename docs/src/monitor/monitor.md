
## Timer

If you want to stop the ODE solver after a certain period of time, you can edit the default option `timer`.

```julia
using RKM
options = SolverOptions(; method = RungeKutta4(), adaptive = Fixed(),
                          timer = TimeLimit(; wtime_minutes = 1),);
```
Here we set the time limit to one minute (or `wtime_minutes = 1`). If you

```julia
function dy_dt!(f, y, t; kwargs...)
    f[1] = (y[1] + 0.5) * (0.5 - y[1])
    sleep(1e-3)
    return nothing
end

t0 = -10.0
y0 = [exp(t0)/(1.0 + exp(t0)) - 0.5]

tf = 10.0
dt0 = 1e-6
```

```julia
julia> @time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options);

┌ Warning: Exceeded time limit of 1.0 minutes (stopping solver...)
└ @ RKM ~/Desktop/RKM.jl/src/timer.jl:119
 60.006475 seconds (104.62 k allocations: 6.267 MiB)
```
*Note: no timer will be set if you used the default value* `wtime_minutes = Inf`.

By default, the solver checks every step whether the runtime has exceeded the limit. You might notice some additional overhead if your ODE function is very fast. If you want to check the runtime less often, you can edit the argument `step_interval`

## Progress bar

```@autodocs
Modules = [RKM]
Pages   = ["src/timer.jl",
           "src/progress.jl"]
```