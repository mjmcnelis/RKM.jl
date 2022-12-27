"""
Curtis' eighth-order method. 

https://link.springer.com/article/10.1007/BF02219778
"""
function Curtis8(; precision::Type{T} = Float64) where T <: AbstractFloat
    # TEMP until can fix get more digits
    if precision == BigFloat || precision == Double64
        @warn "Curtis8 can't use $precision right now (default to Float64)"
        precision = Float64
    end

    butcher = [0 0 0 0 0 0 0 0 0 0 0 0
               1//192 1//192 0 0 0 0 0 0 0 0 0 0
               7.4803184795467909e-2 -4.6236639493684439e-1 5.3716957973231239e-1 0 0 0 0 0 0 0 0 0
               1.1220477719320186e-1 2.8051194298300466e-2 0 8.4153582894901380e-2 0 0 0 0 0 0 0 0
               13//19 8.1897867016409869 0 -2.8775072279672202e1 2.1269496104347006e1 0 0 0 0 0 0 0
               1.7267316464601146e-1 4.0374060815253993e-2 0 0 1.3218823182231848e-1 1.1087200843898359e-4 0 0 0 0 0 0
               8.2732683535398854e-1 2.5873900305712910e-1 0 0 -1.0789704485967246 3.4955683812996113e-1 1.2980014427636228 0 0 0 0 0
               1//2 1.9559695938810262e-1 0 0 -7.7135514072989309e-1 8.8677640841671346e-2 1.0169692666685441 -2.9888726168424973e-2 0 0 0 0
               1.7267316464601146e-1 -5.6453739799316727e-3 0 0 3.6664150851331873e-1 -4.7678144732965888e-2 -2.2046326596199767e-1 1.3899990926683259e-2 6.5918449880904728e-2 0 0 0
               8.2732683535398854e-1 -1.6186994529880696 0 0 8.4859361993123450 -1.6000746209745174 -1.6706963865871344e1 5.6707296863985979e-1 2.6892514721246741 9.0108041351110444 0 0
               1 8.7651711673654731e-1 0 0 -4.1138940838927622 4.7294741782224714e-1 1.3785374631587118e1 2.8529348484093164e-2 -16//27 -9.7062986315222215 2.4941679337757211e-1 0
               1 1//20 0 0 0 0 13//180 1//5 16//45 1//5 13//180 1//20]
    butcher = butcher .|> precision

    RungeKutta(; name = :Curtis_8, butcher)
end

"""
Shanks' eighth-order method. 

https://ntrs.nasa.gov/citations/19650022581
"""
function Shanks8(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0 0 0 0 0 0 0 0 0 0
               1//9 1//9 0 0 0 0 0 0 0 0 0 0 0
               1//6 1//24 3//24 0 0 0 0 0 0 0 0 0 0
               1//4 1//16 0 3//16 0 0 0 0 0 0 0 0 0
               1//10 29//500 0 33//500 -12//500 0 0 0 0 0 0 0 0
               1//6 33//972 0 0 4//972 125//972 0 0 0 0 0 0 0
               1//2 -21//36 0 0 76//36 125//36 -162//36 0 0 0 0 0 0
               2//3 -30//243 0 0 -32//243 125//243 0 99//243 0 0 0 0 0
               1//3 1175//324 0 0 -3456//324 -6250//324 8424//324 242//324 -27//324 0 0 0 0
               5//6 293//324 0 0 -852//324 -1375//324 1836//324 -118//324 162//324 1 0 0 0
               5//6 1303//1620 0 0 -4260//1620 -6875//1620 9990//1620 1030//1620 0 0 1//10 0 0
               1 -8595//4428 0 0 30720//4428 48750//4428 -66096//4428 378//4428 -729//4428 -1944//4428 -1296//4428 3240//4428 0
               1 41//840 0 0 0 0 216//840 272//840 27//840 27//840 36//840 180//840 41//840]
    butcher = butcher .|> precision

    RungeKutta(; name = :Shanks_8, butcher)
end

"""
Shanks' pseudo eighth-order method. 

https://ntrs.nasa.gov/citations/19650022581
"""
function ShanksPseudo8(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0 0 0 0 0 0 0 0
               4//27 4//27 0 0 0 0 0 0 0 0 0
               2//9 1//18 1//6 0 0 0 0 0 0 0 0
               1//3 1//12 0 1//4 0 0 0 0 0 0 0
               1//2 1//8 0 0 3//8 0 0 0 0 0 0
               2//3 13//54 0 -1//2 7//9 4//27 0 0 0 0 0
               1//6 389//4320 0 -1//80 161//720 -103//540 9//160 0 0 0 0
               1 -231//20 0 81//20 -291//5 164//5 -61//10 40 0 0 0
               5//6 -127//288 0 1//16 -113//48 19//12 -1//32 2 1//72 0 0
               1 1481//820 0 -81//820 1776//205 -844//205 18//205 -252//41 -3//41 36//41 0
               1 41//840 0 0 9//280 34//105 9//280 9//35 0 9//35 41//840]
    butcher = butcher .|> precision

    RungeKutta(; name = :Shanks_Pseudo_8, butcher)
end
