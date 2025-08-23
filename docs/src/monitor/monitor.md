
## Timer

You can stop the ODE solver early by changing the default option `timer`

```julia
options = SolverOptions(; method = RungeKutta4(), adaptive = Fixed(),
                          timer = TimeLimit(; wtime_minutes = 1),)
```

where the max number of minutes is set to `wtime_minutes = 1`. 


By default, the solver checks every step whether the runtime has exceeded the limit.

*Note: no timer is set if you used the default value* `wtime_minutes = Inf`.

## Progress bar

```@autodocs
Modules = [RKM]
Pages   = ["src/timer.jl",
           "src/progress.jl"]
```