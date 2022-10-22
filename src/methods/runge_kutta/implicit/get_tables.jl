
# TODO: make docstring
function get_runge_kutta_full_implicit_tables(; precision::Type{<:AbstractFloat})
    (
        LobattoIIIC2(; precision), 
        RadauIA3(; precision), 
        RadauIIA3(; precision),
        
        GaussLegendre4(; precision),
        LobattoIIIA4(; precision), 
        LobattoIIIB4(; precision), 
        LobattoIIIC4(; precision),
        LobattoIIID4(; precision),
    )
end

# TODO: make docstring
function get_runge_kutta_diagonal_implicit_tables(; precision::Type{<:AbstractFloat})
    (
        BackwardEuler1(; precision), 
        ImplicitMidpoint2(; precision), 
        CrankNicolson2(; precision),
        QinZhang2(; precision), 
        KraaijevangerSpijker2(; precision), 
        PareschiRusso2(; precision), 
        LobattoIIIB2(; precision),
        PareschiRusso3(; precision), 
        Crouzeix3(; precision), 
        DIRKL3(; precision),

        Norsett4(; precision),
        LobattoIIICS4(; precision),
    )
end