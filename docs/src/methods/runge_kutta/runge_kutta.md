# Runge-Kutta update

```math 
\frac{dy}{dt} = f(t,y)
```

## Explicit update

```math 
\Delta y_n^{(s)} = \Delta t f\Big(t_n + c_s\Delta t, y_n + \sum_{j=1}^{s-1} A_{sj} \Delta y_n^{(j)}\Big) \\
y_{n+1} = y_n + \sum_{s = 1}^S b_s \Delta y_n^{(s)}
```

```math  
B_{ij} = 1
```