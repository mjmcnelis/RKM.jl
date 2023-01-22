
# Medium-order methods (4-6)

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

- `Fehlberg45()`: Fehlberg's fourth(fifth)-order method
- `CashKarp54()`: Cash and Karp's fifth(fourth)-order method
- `DormandPrince54()`: Dormand and Prince's fifth(fourth)-order method
- `BogackiShampine54()`: Bogacki and Shampine's fifth(fourth)-order method
- `Tsitouras54()`: Bogacki and Shampine's fifth(fourth)-order method
- `Verner56()`: Verner's fifth(sixth)-order method
- `Verner65()`: Verner's sixth(fifth)-order method

*Note: SSP stands for strong stability preserving*

```@autodocs
Modules = [RKM]
Pages   = ["methods/runge_kutta/explicit/fixed/medium_order.jl", 
           "methods/runge_kutta/explicit/embedded/medium_order.jl"]
```