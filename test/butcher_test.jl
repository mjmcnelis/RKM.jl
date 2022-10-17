using Revise
using RKM 
using Test

# TODO: export csv files for larger tables for viewing

butcher_tables = (
    # fixed explicit RK
    Euler1(), Heun2(), Midpoint2(), Ralston2(), Generic2(; alpha = 1), Heun3(), Ralston3(),
    RungeKutta3(), ShuOsher3(),SpiteriRuuth3(), Generic3(; alpha = 1//2),

    RungeKutta4(), ThreeEightsRule4(), Ralston4(), Ketcheson4(), Butcher5(), Butcher6(),

    Curtis8(), Shanks8(), ShanksPseudo8(),

    # embedded explicit RK
    Fehlberg12(), HeunEuler21(), BogackiShampine32(),

    Fehlberg45(), CashKarp54(), DormandPrince54(), BogackiShampine54(), Tsitouras54(), 
    Verner56(), Verner65(),

    Fehlberg78(), DormandPrince87(),

    Feagin108(),

    # fixed implicit RK
    BackwardEuler1(), ImplicitMidpoint2(), CrankNicolson2(), QinZhang2(), 
    KraaijevangerSpijker2(), PareschiRusso2(), LobattoIIIB2(), LobattoIIIC2(), 
    PareschiRusso3(), Crouzeix3(), RadauIA3(), RadauIIA3(), DIRKL3(),

    Norsett4(), LobattoIIIA4(), LobattoIIIB4(), LobattoIIIC4(),
)

debug_table.(butcher_tables)

println("\ndone\n")
