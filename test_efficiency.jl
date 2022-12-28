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

adaptive = Doubling()             
t_span = TimeSpan(; t0 = -10.0, tf = 10.0, dt0 = 1e-4)

t0 = t_span.t0
N = 1
y0 = Float64[]
for i = 1:N
    local a = N == 1 ? 0.5 : 0.5 - 0.25*(i-1.0)/(N-1.0)
    push!(y0, exp(t0) / (1.0 + exp(t0)) - a)
end

methods = OrderedDict(
                    #   Doubling() => [Heun2(), Midpoint2(), Ralston2(), 
                    #                  Heun3(), Ralston3(), RungeKutta3(), ShuOsher3(), 
                    #                  SpiteriRuuth3(), RungeKutta4(), ThreeEightsRule4(),
                    #                  Ralston4(), Ketcheson4(), Butcher5(), Butcher6(),
                    #                  Curtis8(), Shanks8(), ShanksPseudo8()],
                                     
                      Embedded() => [Fehlberg12(), BogackiShampine32(), Fehlberg45(),
                                     CashKarp54(), DormandPrince54(), BogackiShampine54(),
                                     Tsitouras54(), Verner56(), Verner65(), 
                                    #  Fehlberg78(), 
                                     DormandPrince87()
                                     ],
                    )
epsilon_vect = 10.0.^LinRange(-6, -11, 61)

@time plt = efficiency_curve(y0, y_exact, dy_dt!; methods, epsilon_vect, t_span, plot, plot!)

display(plt)

println("\ndone")
