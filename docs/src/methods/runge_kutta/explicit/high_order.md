
# High-order (7-9) explicit Runge-Kutta methods

## Standard

The following standard methods are compatible with `Fixed()`, `Doubling()` and `FiniteDiff()` time-stepping:

- `Curtis8()`: Curtis' eighth-order method
- `Shanks8()`: Shanks' eighth-order method
- `ShanksPseudo8()`: Shanks' pseudo eighth-order method

## Embedded

The following embedded methods are compatible with all adaptive time step options:

- `Fehlberg78()`: Fehlberg's seventh(eighth)-order method
- `DormandPrince87()`: Dormand and Prince's eighth(seventh)-order method

## API Reference

```@autodocs
Modules = [RKM]
Pages   = ["methods/runge_kutta/explicit/fixed/high_order.jl",
           "methods/runge_kutta/explicit/embedded/high_order.jl"]
```