
# TODO: make docstring
function get_runge_kutta_explicit_tables(; precision::Type{T}) where T <: AbstractFloat
    [
        # fixed
        Euler1(; precision),
        Heun2(; precision),
        Midpoint2(; precision),
        Ralston2(; precision),
        Heun3(; precision),
        Ralston3(; precision),
        Kutta3(; precision),
        ShuOsher3(; precision),
        SpiteriRuuth3(; precision),

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
        Fehlberg2(; precision),
        Heun2(; precision),
        BogackiShampine3(; precision),

        Fehlberg5(; precision),
        CashKarp5(; precision),
        DormandPrince5(; precision),
        BogackiShampine5(; precision),
        Tsitouras5(; precision),
        Verner5(; precision),
        Verner6(; precision),

        Fehlberg78(; precision),

        DormandPrince87(; precision),

        Feagin108(; precision),

        Feagin1210(; precision),

        Feagin1412(; precision),
    ]
end
