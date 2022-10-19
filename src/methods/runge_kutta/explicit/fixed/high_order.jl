
function Curtis8(; precision::Type{<:AbstractFloat} = Float64)
    # TEMP until can fix get more digits
    if precision == BigFloat
        @warn "Curtis8 can't use BigFloat right now (default to Float64)"
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

function ShanksPseudo8(; precision::Type{<:AbstractFloat} = Float64)
    # TEMP until can fix get more digits
    if precision == BigFloat
        @warn "ShanksPseudo8 can't use BigFloat right now (default to Float64)"
        precision = Float64
    end

    # TODO: get fractions
    butcher = [0 0 0 0 0 0 0 0 0 0 0
               0.14814814814814814 0.14814814814814814 0 0 0 0 0 0 0 0 0
               2//9 1//18 1//6 0 0 0 0 0 0 0 0
               1//3 1//12 0 1//4 0 0 0 0 0 0 0
               1//2 1//8 0 0 3//8 0 0 0 0 0 0
               0.6666666666666666 0.24074074074074073 0 -0.5 0.7777777777777778 0.14814814814814814 0 0 0 0 0
               0.16666666666666666 0.09004629629629629 0 -0.0125 0.22361111111111112 -0.19074074074074074 0.05625 0 0 0 0
               1 -11.55 0 4.05 -58.2 32.8 -6.1 40 0 0 0
               0.8333333333333334 -0.4409722222222222 0 0.0625 -2.3541666666666665 1.5833333333333333 -0.03125 2 0.013888888888888888 0 0
               1 1.8060975609756098 0 -0.09878048780487805 8.663414634146342 -4.117073170731707 0.08780487804878048 -6.146341463414634 -0.07317073170731707 0.8780487804878049 0
               1 0.04880952380952381 0 0 0.03214285714285714 0.3238095238095238 0.03214285714285714 0.2571428571428571 0 0.2571428571428571  0.04880952380952381]
    butcher = butcher .|> precision

    RungeKutta(; name = :Shanks_Pseudo_8, butcher)
end

    