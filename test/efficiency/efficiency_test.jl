using Revise
using RKM
import DoubleFloats: Double64
import DataStructures: OrderedDict
import Plots: plotly
plotly()
include("$RKM_root/test/efficiency/efficiency_curve.jl")
!(@isdefined dy_dt!) ? include("$RKM_root/equations.jl") : nothing 

precision = Double64 
# precision = BigFloat

t_range = TimeRange(; t0 = -10.0, tf = 10.0, dt0 = 1e-4)

const C = 0.5 
t0 = t_range.t0 |> precision
y0 = exp(t0)/(1.0 + exp(t0)) - C 

methods = OrderedDict(
    Doubling() => [ 
                   Heun2(), Midpoint2(), Ralston2(), 
                   Heun3(), Ralston3(), RungeKutta3(), ShuOsher3(), SpiteriRuuth3(), 
                   RungeKutta4(), ThreeEightsRule4(), Ralston4(), Ketcheson4(), 
                   Butcher5(), 
                   Butcher6(),
                   Curtis8(), Shanks8(), ShanksPseudo8(),
                ],      
    Embedded() => [
                   HeunEuler21(),
                   BogackiShampine32(), 
                   Fehlberg45(),
                   CashKarp54(), DormandPrince54(), BogackiShampine54(), Tsitouras54(), Verner56(),
                   Verner65(), 
                   Fehlberg78(), 
                   DormandPrince87(),
                   Feagin108(),
                   Feagin1210(),
                   Feagin1412(),
                ],
)
epsilon_vect = 10.0.^LinRange(-4, -28, 25)
# farthest I got with Feagin1412
# epsilon_vect = 10.0.^LinRange(-4, -80, 77)

@time plt = efficiency_curve(y0, y_exact, dy_dt!; precision, methods, 
                             epsilon_vect, t_range)

display(plt)

println("\ndone")
