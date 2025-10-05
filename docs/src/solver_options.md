
# Solver options

The ODE evolution function `evolve_ode` requires a set of solver options, which are stored in the struct `SolverOptions`. Here is an example set of all the available options:

```julia
options = SolverOptions(;
    method = RungeKutta4(),
    adaptive = Fixed(),
    timer = TimeLimit(),
    state_jacobian = FiniteJacobian(),
    root_finder = Newton(),
    eigenmax = NoEigenMax(),
    sensitivity = NoSensitivity(),
    interpolator = NoInterpolation(),
    save_solution = true
    save_time_derivative = false
    show_progress = false
    time_subroutine = false
    precision = Float64
)
```

We provide an overview of each solver option:

##### 1. `method`
You must specify the ODE method used for the time evolution. The solver supports two classes of ODE methods: Runge--Kutta and linear multistep.
    - `Euler1()`
    - sgdgs
```julia
julia> list_explicit_runge_kutta_methods()
Low order (1-3)       | Euler1, Heun2, Midpoint2, Ralston2, Fehlberg2, Heun2, Heun3,
                      | Ralston3, Kutta3, ShuOsher3, SpiteriRuuth3, BogackiShampine3
-------------------------------------------------------------------------------
Medium order (4-6)    | RungeKutta4, ThreeEightsRule4, Ralston4, Ketcheson4, Butcher5,
                      | Fehlberg5, CashKarp5, DormandPrince5, BogackiShampine5,
                      | Tsitouras5, Verner5, Butcher6, Verner6
-------------------------------------------------------------------------------
High order (7-9)      | Fehlberg7, DormandPrince8, Curtis8, Shanks8, ShanksPseudo8
-------------------------------------------------------------------------------
Very high order (10+) | Feagin10, Feagin12, Feagin14

julia> list_diagonal_implicit_runge_kutta_methods()
Low order (1-3)       | BackwardEuler1, TrapezoidRuleBDF2, ImplicitTrapezoid2,
                      | ImplicitMidpoint2, QinZhang2, KraaijevangerSpijker2,
                      | PareschiRusso2, LobattoIIIB2, PareschiRusso3, Crouzeix3,
                      | DIRKL3
-------------------------------------------------------------------------------
Medium order (4-6)    | Norsett4, LobattoIIICS42
```
##### 2. `adaptive`
...

```@autodocs
Modules = [RKM]
Pages   = ["src/options.jl"]
```