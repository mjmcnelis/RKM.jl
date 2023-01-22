
# Method construction

In this section, we discuss how a Runge-Kutta method's Butcher tableau and properties
are constructed and stored in the `RungeKutta` struct.

## Example 

Here we set the ODE solver method to the classic fourth-order Runge-Kutta method:
```julia
method = RungeKutta4()
```

The code for this method's constructor is

```julia 
function RungeKutta4(; precision::Type{T} = Float64) where T <: AbstractFloat
    butcher = [0 0 0 0 0
               1//2 1//2 0 0 0
               1//2 0 1//2 0 0
               1 0 0 1 0
               1 1//6 1//3 1//3 1//6]
    butcher = butcher .|> precision

    return RungeKutta(; name = :Runge_Kutta_4, butcher)
end
```
First we write out the Butcher tableau in `Matrix` form. If the coefficients are simple,
we express them as integers and fractions (otherwise we use `big` strings). The matrix is 
converted to the type `precision`, which is set by the user. 

*Note: the `1` entry in the lower-left corner of the Butcher tableau is only used to verify that the `b` coefficients sum up to unity within floating precision.*
## API Reference 

```@autodocs
Modules = [RKM]
Pages   = ["methods/runge_kutta/runge_kutta.jl", 
           "methods/properties.jl", 
           "methods/code_names.jl"]
```