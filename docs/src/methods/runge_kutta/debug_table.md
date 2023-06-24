
# Debugging the Butcher tableau

## Testing the order conditions

The test file `test/butcher_test.jl` can be used to debug the Butcher tableau of each
Runge-Kutta method to a certain extent. It checks whether or not the order conditions are
satisfied within numerical precision:
```math
c_i = \sum_{j=1}^{S} A_{ij} \,\, \forall \,\, i \\
\sum_{j=1}^{S} b_j = 1 \\
\sum_{j=1}^{S} \hat{b}_j = 1 \,,
```
where `S` are the number of stages. Any individual stage that fails the test will register
as either `Broken` or `Fail`, depending on how badly the condition is violated. This test
is particularly important when using high-order methods with many coefficients and a high
numerical precision is needed. To demonstrate, we validate the order conditions of the
Feagin 14(12) method up to double-float precision:
```julia
julia> using RKM
julia> import DoubleFloats: Double64
julia> method = Feagin1412(; precision = Double64);
julia> debug_table(method)
#┌ Warning: |B[18,1] - ∑_{j>1} B[18,j]| = 7.7e-34 > 4.84e-34. Check row 18 in Feagin1412.
#└ @ RKM ~/Desktop/RKM.jl/src/methods/runge_kutta/debug_table.jl:50
#┌ Warning: |B[20,1] - ∑_{j>1} B[20,j]| = 2.1e-31 > 1.8599999999999997e-31. Check row 20 in Feagin1412.
#└ @ RKM ~/Desktop/RKM.jl/src/methods/runge_kutta/debug_table.jl:50
```
We see that two stages are broken, . However. We caution that this is.
It may not detect a faulty row where the errors of
multiple coefficients cancel out. defects in multiple coefficients per
row whose errors cancel out.


## API Reference

```@autodocs
Modules = [RKM]
Pages   = ["methods/runge_kutta/debug_table.jl"]
```