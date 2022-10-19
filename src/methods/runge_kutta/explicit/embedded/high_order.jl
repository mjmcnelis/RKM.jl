
function Fehlberg78(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0 0 0 0 0 0 0 0 0 0 0
               2//27 2//27 0 0 0 0 0 0 0 0 0 0 0 0
               1//9 1//36 1//12 0 0 0 0 0 0 0 0 0 0 0
               1//6 1//24 0 1//8 0 0 0 0 0 0 0 0 0 0
               5//12 5//12 0 -25//16 25//16 0 0 0 0 0 0 0 0 0
               1//2 1//20 0 0 1//4 1//5 0 0 0 0 0 0 0 0
               5//6 -25//108 0 0 125//108 -65//27 125//54 0 0 0 0 0 0 0
               1//6 31//300 0 0 0 61//225 -2//9 13//900 0 0 0 0 0 0
               2//3 2 0 0 -53//6 704//45 -107//9 67//90 3 0 0 0 0 0
               1//3 -91//108 0 0 23//108 -976//135 311//54 -19//60 17//6 -1//12 0 0 0 0
               1 2383//4100 0 0 -341//164 4496//1025 -301//82 2133//4100 45//82 45//164 18//41 0 0 0
               0 3//205 0 0 0 0 -6//41 -3//205 -3//41 3//41 6//41 0 0 0
               1 -1777//4100 0 0 -341//164 4496//1025 -289//82 2193//4100 51//82 33//164 12//41 0 1 0
               1 41//840 0 0 0 0 34//105 9//35 9//35 9//280 9//280 41//840 0 0
               1 0 0 0 0 0 34//105 9//35 9//35 9//280 9//280 0 41//840 41//840]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Fehlberg_78, butcher)
end

function DormandPrince87(; precision::Type{<:AbstractFloat} = Float64)
    # TEMP until can fix large Ints
    if precision == BigFloat
        @warn "DormandPrince87 can't use BigFloat right now (default to Float64)"
        precision = Float64
    end

    butcher = [0 0 0 0 0 0 0 0 0 0 0 0 0 0
               1//18 1//18 0 0 0 0 0 0 0 0 0 0 0 0
               1//12 1//48 1//16 0 0 0 0 0 0 0 0 0 0 0
               1//8 1//32 0 3//32 0 0 0 0 0 0 0 0 0 0
               5//16 5//16 0 -75//64 75//64 0 0 0 0 0 0 0 0 0
               3//8 3//80 0 0 3//16 3//20 0 0 0 0 0 0 0 0
               59//400 29443841//614563906 0 0 77736538//692538347 -28693883//1125000000 23124283//1800000000 0 0 0 0 0 0 0
               93//200 16016141//946692911 0 0 61564180//158732637 22789713//633445777 545815736//2771057229 -180193667//1043307555 0 0 0 0 0 0
               5490023248//9719169821 39632708//573591083 0 0 -433636366//683701615 -421739975//2616292301 100302831//723423059 790204164//839813087 800635310//3783071287 0 0 0 0 0
               13//20 246121993//1340847787 0 0 -37695042795//15268766246 -309121744//1061227803 -12992083//490766935 6005943493//2108947869 393006217//1396673457 123872331//1001029789 0 0 0 0
               1201146811//1299019798 -1028468189//846180014 0 0 8478235783//508512852 1311729495//1432422823 -10304129995//1701304382 -48777925059//3047939560 15336726248//1032824649 -45442868181//3398467696 3065993473//597172653 0 0 0
               1 185892177//718116043 0 0 -3185094517//667107341 -477755414//1098053517 -703635378//230739211 5731566787//1027545527 5232866602//850066563 -4093664535//808688257 3962137247//1805957418 65686358//487910083 0 0
               1 403863854//491063109 0 0 -5068492393//434740067 -411421997//543043805 652783627//914296604 11173962825//925320556 -13158990841//6184727034 3936647629//1978049680 -160528059//685178525 248638103//1413531060 0 0
               1 13451932//455176623 0 0 0 0 -808719846//976000145 1757004468//5645159321 656045339//265891186 -3867574721//1518517206 465885868//322736535  53011238//667516719 2//45 0
               1 14005451//335480064 0 0 0 0 -59238493//1068277825 181606767//758867731 561292985//797845732 -1041891430//1371343529 760417239//1151165299 118820643//751138087 -528747749//2220607170 1//4]
#    @show typeof(butcher)
#    q()
    butcher = butcher .|> precision 

    RungeKutta(; name = :Dormand_Prince_87, butcher)
end
