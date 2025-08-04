
# Low-order (1-3) implicit Runge-Kutta methods

## Standard

Currently, only the standard diagonal implicit methods are compatible with `Fixed()`
time-stepping.

- `BackwardEuler1()`: first-order backward Euler method (L-stable)
- `ImplicitTrapezoid2()`: second-order implicit trapezoid rule (A-stable)
- `ImplicitMidpoint2()`: second-order implicit mid-point rule (symplectic, A-stable)
- `QinZhang2()`: Qin and Zhang's second-order method
- `KraaijevangerSpijker2()`: Kraaijevanger and Spijker's second-order method
- `PareschiRusso2()`: Pareschi and Russo's second-order method
- `PareschiRusso3()`: Pareschi and Russo's third-order method
- `Crouzeix3()`: Crouzeix's third-order method
- `RadauIA3()`: Radau IA3 third-order method
- `RadauIIA3()`: Radau IA3 third-order method
- `DIRKL3()`: a third-order L-stable diagonal implicit method (name unknown)

## Embedded

The following embedded methods are compatible with all adaptive time step options:

- `Fehlberg2()`: Fehlberg's second-order method
- `Heun2()`: Heun-Euler second-order method
- `BogackiShampine3()`: Bogacki and Shampine's third-order method

*Note: SSP stands for strong stability preserving*

## API Reference

```@autodocs
Modules = [RKM]
Pages   = ["methods/runge_kutta/implicit/fixed/low_order.jl",
           "methods/runge_kutta/implicit/embedded/low_order.jl"]
```