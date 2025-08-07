
# Medium-order (4-6) explicit Runge-Kutta methods

## Standard

The following standard methods are compatible with `Fixed()`, `Doubling()` and `FiniteDiff()` time-stepping:

- `RungeKutta4()`: Runge and Kutta's classic fourth-order method
- `ThreeEightsRule4()`: fourth-order 3/8 rule
- `Ralston4()`: Ralston's fourth-order method
- `Ketcheson4()`: Ketcheson's fourth-order SSP method
- `Butcher5()`: Butcher's fifth-order method
- `Butcher6()`: Butcher's sixth-order method

## Embedded

The following embedded methods are compatible with all adaptive time step options:

- `Fehlberg5()`: Fehlberg's fifth-order method
- `CashKarp5()`: Cash and Karp's fifth-order method
- `DormandPrince5()`: Dormand and Prince's fifth-order method
- `BogackiShampine5()`: Bogacki and Shampine's fifth-order method
- `Tsitouras5()`: Bogacki and Shampine's fifth-order method
- `Verner5()`: Verner's fifth-order method
- `Verner6()`: Verner's sixth-order method

*Note: SSP stands for strong stability preserving*

## API Reference

```@autodocs
Modules = [RKM]
Pages   = ["methods/runge_kutta/explicit/fixed/medium_order.jl",
           "methods/runge_kutta/explicit/embedded/medium_order.jl"]
```