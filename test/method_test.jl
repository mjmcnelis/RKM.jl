using Revise
using RKM 

# TODO: make butcher tableau debugger

fixed_explicit = (
    Euler1(), Heun2(), Midpoint2(), Ralston2(),
    Heun3(), Ralston3(), RungeKutta3(), ShuOsher3(), SpiteriRuuth3(),
    RungeKutta4(), ThreeEightsRule4(), Ralston4(), Ketcheson4(), 
    Butcher5(), Butcher6(),
    Curtis8(), Shanks8(), ShanksPseudo8()
)

embedded_explicit =  (
    Fehlberg12(), HeunEuler21(), BogackiShampine32()
)


# TODO: export csv files for larger tables
for x in fixed_explicit println(x) end
for x in embedded_explicit println(x) end
