
# TODO: make docstring
function get_runge_kutta_explicit_tables(; precision::Type{T}) where T <: AbstractFloat
    [
        # fixed
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

        # embedded
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
        
        Feagin1210(; precision),
    ]
end
