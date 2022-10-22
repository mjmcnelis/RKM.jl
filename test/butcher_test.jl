using Revise
using RKM 
using Test
# TODO: export csv files for larger tables for viewing

function get_butcher_tables(; precision::Type{<:AbstractFloat})
    (
        # fixed explicit RK
        Euler1(; precision), 
        Heun2(; precision), 
        Midpoint2(; precision), 
        Ralston2(; precision), 
        Generic2(; alpha = 1, precision), 
        Heun3(; precision),
        Ralston3(; precision),
        RungeKutta3(; precision), 
        ShuOsher3(; precision), 
        SpiteriRuuth3(; precision), 
        Generic3(; alpha = 1//2, precision),

        RungeKutta4(; precision), 
        ThreeEightsRule4(; precision), 
        Ralston4(; precision), 
        Ketcheson4(; precision), 
        Butcher5(; precision), 
        Butcher6(; precision),

        Curtis8(; precision), 
        Shanks8(; precision), 
        ShanksPseudo8(; precision),

        # embedded explicit RK
        Fehlberg12(; precision), 
        HeunEuler21(; precision), 
        BogackiShampine32(; precision),

        Fehlberg45(; precision), 
        CashKarp54(; precision),
        DormandPrince54(; precision), 
        BogackiShampine54(; precision), 
        Tsitouras54(; precision), 
        Verner56(; precision), 
        Verner65(; precision),

        Fehlberg78(; precision), 
        DormandPrince87(; precision),

        Feagin108(; precision),

        # fixed implicit RK
        BackwardEuler1(; precision), 
        ImplicitMidpoint2(; precision), 
        CrankNicolson2(; precision), 
        QinZhang2(; precision), 
        KraaijevangerSpijker2(; precision), 
        PareschiRusso2(; precision), 
        LobattoIIIB2(; precision), 
        LobattoIIIC2(; precision), 
        PareschiRusso3(; precision), 
        Crouzeix3(; precision), 
        RadauIA3(; precision), 
        RadauIIA3(; precision), 
        DIRKL3(; precision),

        Norsett4(; precision), 
        LobattoIIIA4(; precision), 
        LobattoIIIB4(; precision), 
        LobattoIIIC4(; precision)
    )
end

# butcher_tables = get_butcher_tables(; precision = Float32)
# butcher_tables = get_butcher_tables(; precision = Float64)
butcher_tables = get_butcher_tables(; precision = BigFloat)

debug_table.(butcher_tables)

println("\ndone\n")
