
# Low-order (1-3) explicit Runge-Kutta methods

## Standard

The following standard methods are compatible with `Fixed()`, `Doubling()` and `FiniteDiff()` time-stepping:

- `Euler1()`: Euler's first-order method
- `Heun2()`: Heun's second-order SSP method
- `Midpoint2()`: second-order mid-point rule
- `Ralston2()`: Ralston's second-order method
- `Heun3()`: Heun's third-order method
- `Ralston3()`: Ralston's third-order
- `Kutta3()`: Kutta's third-order method
- `ShuOsher3()`: Shu and Osher's third-order SSP method
- `SpiteriRuuth3()`: Spiteri and Ruuth's third-order SSP method

## Embedded

The following embedded methods are compatible with all adaptive time step options:

- `Fehlberg21()`: Fehlberg's second(first)-order method
- `HeunEuler21()`: Heun-Euler second(first)-order method
- `BogackiShampine32()`: Bogacki and Shampine's third(second)-order method

*Note: SSP stands for strong stability preserving*

## API Reference

```@autodocs
Modules = [RKM]
Pages   = ["methods/runge_kutta/explicit/fixed/low_order.jl",
           "methods/runge_kutta/explicit/embedded/low_order.jl"]
```