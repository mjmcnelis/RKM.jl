
# Solver options

The ODE evolution function `evolve_ode` requires a set of solver options, which are stored in the struct `SolverOptions`. Here is an example set of all the options:

```julia
options = SolverOptions(;
    method = RungeKutta4(),
    adaptive = Fixed(),
    timer = TimeLimit(),
    state_jacobian = FiniteJacobian(),
    root_finder = Newton(),
    eigenmax = NoEigenMax(),
    sensitivity = NoSensitivity(),
    interpolator = NoInterpolation(),
    save_solution = true
    save_time_derivative = false
    show_progress = false
    time_subroutine = false
    precision = Float64
)
```

We provide an overview of each solver option:

## `method`

The ODE method used for the time evolution must be specified.

The solver supports two classes of ODE methods: Runge--Kutta and linear multistep. The available Runge--Kutta methods are either explicit or diagonal-implicit (full-implicit methods are not supported yet). Two examples include the classic RK4 method (explicit)

```julia
method = RungeKutta4()
```

and the TRBDF2 method (diagonal-implicit), which is robust for stiff ODE problems

```julia
method = TrapezoidRuleBDF2()
```

A full list of Runge--Kutta methods can be found in the pages [Explicit Runge-Kutta methods](methods/runge_kutta/explicit_runge_kutta.md) and [Implicit Runge-Kutta methods](methods/runge_kutta/implicit_runge_kutta.md).

*Note: you do not have to pass a separate* `precision` *argument to the method's outer constructor (e.g.* `RungeKutta4(; precision = Float64)`*); the solver reconstructs the ODE method with the*
`precision` *field set in the solver options.*

!!! warning "TODO: mention multistep after reorg (e.g. AdamsBashforth2) and more support"

## `adaptive`

The method used to compute the adaptive time step must be specified.

If you want the time step to be constant, set

```julia
adaptive = Fixed()
```

This fixes the time step to the initial value `dt0` passed to the `evolve_ode` function. Otherwise if you want an adaptive time step, you can choose either step doubling and embedded.

```julia
adaptive = Doubling()   # or Embedded()
```

All Runge--Kutta methods are compatible with step doubling, but it is only practical for methods with a low number of stages (1 or 2).

## `save_solution`

The numerical ODE solution is stored if set to `true` (defaulted to `true`).

In practice, you will want to save the ODE solution, but setting `save_solution = false` can be useful in dry run testing.

## `save_time_derivative`

The first-order time derivatives are also stored if set to `true` (defaulted to `false`).

See the [Time derivatives](solution/solution_data.md#Time-derivatives) section for an example.

*Note: no time derivative data can be saved if `save_solution = false`.*

## `show_progress`

A progress bar is displayed if set to `true` (defaulted to `false`).

For more details, go to the [Progress bar](monitor/monitor.md#Progress-bar) section.

## `time_subroutine`

Core subroutines in the time evolution loop are timed if set to `true` (defaulted to `false`).

For more details, go to the [Subroutine times](statistics.md#Subroutine-times) section.

## `precision`

The floating point type used in the solver (defaulted to `Float64`).

You can increase the floating precision for ODE problems that require high numerical accuracy. The following arbitrary precision types are supported:

- `Double64` from `DoubleFloats.jl` (31 digits)
- `BigFloat` (76 digits, default)

The float type used for the initial conditions `(t0, y0)` is independent of `precision`.

For implicit ODE methods, we recommend using the default value `precision = Float64`.

## In-place evolution

```julia
sol = Solutirecision` to the on(options)
evolve_ode!(sol, ...)
```
## API reference

```@autodocs
Modules = [RKM]
Pages   = ["src/options.jl"]
```