using Revise
using RKM 
using Test

# TODO: export csv files for larger tables for viewing 

butcher_tables = (
    # fixed explicit
    Euler1(), Heun2(), Midpoint2(), Ralston2(), Generic2(; alpha = 1),
    Heun3(), Ralston3(), RungeKutta3(), ShuOsher3(), 
    SpiteriRuuth3(), Generic3(; alpha = 1//2),
    RungeKutta4(), ThreeEightsRule4(), Ralston4(), Ketcheson4(), 
    Butcher5(), Butcher6(),
    Curtis8(), Shanks8(), ShanksPseudo8(),
    # embedded explicit
    Fehlberg12(), HeunEuler21(), BogackiShampine32(),
    Fehlberg45(), CashKarp54(), DormandPrince54(), BogackiShampine54(),
    Tsitouras54(), Verner56(), Verner65(),
    Fehlberg78(), DormandPrince87(),
    Feagin108(),
)

for x in butcher_tables
    @test debug_table(x) == true
end 

println("\ndone")
