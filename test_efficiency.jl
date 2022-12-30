using Revise
using RKM
import DoubleFloats: Double64
import DataStructures: OrderedDict
import Plots: plot, plot!, plotly
plotly()
try
    dy_dt!
    y_exact
catch err
    isa(err, UndefVarError) ? include("$RKM_root/equations.jl") : nothing
end
            
t_span = TimeSpan(; t0 = -10.0, tf = 10.0, dt0 = 1e-4)

# precision = Double64 
precision = BigFloat 

t0 = t_span.t0 |> precision
N = 1
y0 = precision[]
for i = 1:N
    local a = N == 1 ? 0.5 : 0.5 - 0.25*(i-1.0)/(N-1.0)
    push!(y0, exp(t0) / (1.0 + exp(t0)) - a)
end

methods = OrderedDict(
    Doubling() => [Heun2(; precision), Midpoint2(; precision), Ralston2(; precision), 
                   Heun3(; precision), Ralston3(; precision), RungeKutta3(; precision),
                   ShuOsher3(; precision), SpiteriRuuth3(; precision), 
                   RungeKutta4(; precision), ThreeEightsRule4(; precision), 
                   Ralston4(; precision), Ketcheson4(; precision), 
                   Butcher5(; precision), 
                   Butcher6(; precision),
                   Curtis8(; precision), Shanks8(; precision), ShanksPseudo8(; precision),
                ],      
    Embedded() => [HeunEuler21(; precision),
                   BogackiShampine32(; precision), 
                   Fehlberg45(; precision),
                   CashKarp54(; precision), DormandPrince54(; precision), 
                   BogackiShampine54(; precision), Tsitouras54(; precision), 
                   Verner56(; precision), 
                   Verner65(; precision), 
                   DormandPrince87(; precision),
                   Feagin108(; precision),
                ],
)
epsilon_vect = 10.0.^LinRange(-4, -17, 14)

@time plt = efficiency_curve(y0, y_exact, dy_dt!; methods, epsilon_vect, t_span, plot, plot!)

display(plt)

println("\ndone")
