using Revise
using RKM 

# TODO: make butcher tableau debugger

RK_methods = (
    Euler1(), Heun2(), Midpoint2(), Ralston2(),
    Heun3(), Ralston3(), RungeKutta3(), ShuOsher3(), SpiteriRuuth3(),
    RungeKutta4(), ThreeEightsRule4(), Ralston4(), Ketcheson4()
)

for x in RK_methods println(x) end
