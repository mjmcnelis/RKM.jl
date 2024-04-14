
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
function RungeKutta4(precision::Type{T} = Float64) where T <: AbstractFloat
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
we express them as integers and fractions (otherwise we use `big` strings). The matrix
`butcher` is converted to the float type `precision`, which is set by the user.

Next, we pass `butcher` and a symbolic label `name` to `RungeKutta`'s outer constructor.
The outer constructor partitions `butcher` into the arrays `c`, `A_T` and `b` (`b_hat`)
and stores them as static types. It also infers several properties, such as the number of
stages and whether the method is explicit or implicit. The order of the primary (embedded)
update is extracted by parsing `name`. The code below shows the fields of the `RungeKutta`
struct `method` after it is set to `RungeKutta4()`:

```julia
julia> for field in method |> typeof |> fieldnames
           println("$field = $(getproperty(method, field))")
       end
name = Runge_Kutta_4
c = [0.0, 0.5, 0.5, 1.0]
A_T = [0.0 0.5 0.0 0.0; 0.0 0.0 0.5 0.0; 0.0 0.0 0.0 1.0; 0.0 0.0 0.0 0.0]
b = [0.16666666666666666, 0.3333333333333333, 0.3333333333333333, 0.16666666666666666]
b_hat = [0.16666666666666666, 0.3333333333333333, 0.3333333333333333, 0.16666666666666666]
stages = 4
order = [4.0]
iteration = Explicit()
fsal = false
code_name = RK4
```

*Note: the `1` entry in the lower-left corner of `butcher` is only used as a placeholder to verify that the coefficients in `b` (`b_hat`) sum up to unity within floating precision.*

## API Reference

```@autodocs
Modules = [RKM]
Pages   = ["methods/runge_kutta/runge_kutta.jl",
           "methods/properties.jl",
           "methods/code_names.jl"]
```