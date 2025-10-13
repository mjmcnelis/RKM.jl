
## Implicit Runge-Kutta methods

### Low-order (1-3)

- `BackwardEuler1()`
- `ImplicitTrapezoid2()`
- `ImplicitMidpoint2()`
- `QinZhang2()`
- `KraaijevangerSpijker2()`
- `PareschiRusso2()`
- `LobattoIIIB2()`
- `PareschiRusso3()`
- `Crouzeix3()`
- `DIRKL3()`

### Medium order (4-6)

- `Norsett4()`

## API Reference

```@autodocs
Modules = [RKM]
Pages   = ["methods/runge_kutta/implicit/low_order.jl",
           "methods/runge_kutta/implicit/fixed/medium_order.jl",
           "methods/runge_kutta/implicit/embedded/medium_order.jl"]
```