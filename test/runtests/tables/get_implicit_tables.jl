

# TODO: make docstring
function get_runge_kutta_diagonal_implicit_tables(; precision::Type{<:AbstractFloat})
    [
        # fixed
        BackwardEuler1(; precision),
        ImplicitMidpoint2(; precision),
        QinZhang2(; precision),
        KraaijevangerSpijker2(; precision),
        PareschiRusso2(; precision),
        PareschiRusso3(; precision),
        Crouzeix3(; precision),
        DIRKL3(; precision),

        Norsett4(; precision),

        # embedded
        CrankNicolson21(; precision),
        LobattoIIIB21(; precision),

        LobattoIIICS42(; precision),
    ]
end

# TODO: make docstring
function get_runge_kutta_full_implicit_tables(; precision::Type{<:AbstractFloat})
    [
        # fixed
        RadauIA3(; precision),
        RadauIIA3(; precision),
        RaduaIA5(; precision),

        #embedded
        LobattoIIIC21(; precision),

        GaussLegendre42(; precision),
        LobattoIIIA42(; precision),
        LobattoIIIB42(; precision),
        LobattoIIIC42(; precision),
        LobattoIIID42(; precision),
        RaduaIIA52(; precision),
        GaussLegendre64(; precision),
    ]
end
