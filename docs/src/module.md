# Module

## Exports 
- `evolve_ode` for solving the ODE system and outputting the solution
- `get_solution` for post-processing the ODE solution
- `Parameters` for storing the ODE solver options
- `TimeRange` for specifying the time evolution interval
- `TimeLimit` for setting a timer for the solver routine
- Runge--Kutta methods (see [Runge--Kutta](methods/runge_kutta/runge_kutta.html))
- Adaptive time step options (see ...make a page...)

## Internal modules
- `StaticArrays`: `SVector`, `SMatrix`, `MVector` and `MMatrix` for static allocations
- `FastBroadcast`: `@..` for fast broadcasting of recurring element-wise operations
- `ProgressMeter`: `Progress` and `next!` for monitoring the status of solver routines

## External modules
- `DoubleFloat`: `Double64` for runs that require double-float precision
