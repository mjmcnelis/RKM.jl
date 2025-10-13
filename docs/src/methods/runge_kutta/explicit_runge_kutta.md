
## Explicit Runge-Kutta methods

### Low-order (1-3)
- `Euler1()`
- `Heun2()`
- `Midpoint2()`
- `Ralston2()`
- `Fehlberg2()`
- `Heun3()`
- `Ralston3()`
- `Kutta3()`
- `ShuOsher3()`
- `SpiteriRuuth3()`
- `BogackiShampine3()`

### Medium-order (4-6)
- `RungeKutta4()`
- `ThreeEightsRule4()`
- `Ralston4()`
- `Ketcheson4()`
- `Butcher5()`
- `Fehlberg5()`
- `CashKarp5()`
- `DormandPrince5()`
- `BogackiShampine5()`
- `Tsitouras5()`
- `Verner5()`
- `Butcher6()`
- `Verner6()`
- `Butcher6()`

### High order (7-9)
- `Fehlberg7()`
- `DormandPrince8()`
- `Curtis8()`
- `Shanks8()`
- `ShanksPseudo8()`

### Very high order (10+)
- `Feagin10()`
- `Feagin12()`
- `Feagin14()`

## API Reference

```@autodocs
Modules = [RKM]
Pages   = ["methods/runge_kutta/explicit/low_order.jl",
           "methods/runge_kutta/explicit/medium_order.jl",
           "methods/runge_kutta/explicit/high_order.jl",
           "methods/runge_kutta/explicit/very_high_order.jl"]
```