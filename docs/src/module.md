# Module

## Exports 
- `evolve_ode` for solving the ODE system and outputting the solution
- `get_solution` for post-processing the ODE solution
- `Parameters` for storing the ODE solver options
- `TimeRange` for specifying the time evolution interval
- `TimeLimit` for setting a timer for the solver routine
- Runge--Kutta methods (see [Runge--Kutta](methods/runge_kutta/runge_kutta.html))
- Adaptive time step options (see ...make a page...)

## Dependencies
- [`StaticArrays.jl`](https://github.com/JuliaArrays/StaticArrays.jl): `SVector`, `SMatrix`, `MVector`, `MMatrix`, `@MVector` and `@MMatrix` for the static allocation of Butcher tableaus, state vectors and intermediate caches
- [`FastBroadcast.jl`](https://github.com/YingboMa/FastBroadcast.jl): `@..` for broadcasting recurring element-wise operations in the state update routines
- [`ProgressMeter.jl`](https://github.com/timholy/ProgressMeter.jl): `Progress` and `next!` for monitoring the real-time status of the solver routine

## External modules
- [`DoubleFloats.jl`](https://github.com/JuliaMath/DoubleFloats.jl): `Double64` for runs that require quadruple-float precision
- [`BenchmarkTools.jl`](https://github.com/JuliaCI/BenchmarkTools.jl): `@benchmark`, `@btime` and `mean` for benchmarking solver runtime and memory usage
