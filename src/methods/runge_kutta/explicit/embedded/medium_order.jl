
function Fehlberg45(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0 0 0 0
               1//4 1//4 0 0 0 0 0
               3//8 3//32 9//32 0 0 0 0
               12//13 1932//2197 -7200//2197 7296//2197 0 0 0
               1 439//216 -8 3680//513 -845//4104 0 0
               1//2 -8//27 2 -3544//2565 1859//4104 -11//40 0
               1 25//216 0 1408//2565 2197//4104 -1//5 0
               1 16//135 0 6656//12825 28561//56430 -9//50 2//55]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Fehlberg_45, butcher)
end

function CashKarp54(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0 0 0 0
               1//5 1//5 0 0 0 0 0
               3//10 3//40 9//40 0 0 0 0
               3//5 3//10 -9//10 6//5 0 0 0
               1 -11//54 5//2 -70//27 35//27 0 0
               7//8 1631//55296 175//512 575//13824 44275//110592 253//4096 0
               1 37//378 0 250//621 125//594 0 512//1771
               1 2825//27648 0 18575//48384 13525//55296 277//14336 1//4]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Cash_Karp_54, butcher)
end

function DormandPrince54(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0 0 0 0 0
               1//5 1//5 0 0 0 0 0 0
               3//10 3//40 9//40 0 0 0 0 0
               4//5 44//45 -56//15 32//9 0 0 0 0
               8//9 19372//6561 -25360//2187 64448//6561 -212//729 0 0 0
               1 9017//3168 -355//33 46732//5247 49//176 -5103//18656 0 0
               1 35//384 0 500//1113 125//192 -2187//6784 11//84 0
               1 35//384 0 500//1113 125//192 -2187//6784 11//84 0
               1 5179//57600 0 7571//16695 393//640 -92097//339200 187//2100 1//40]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Dormand_Prince_54, butcher)
end

function BogackiShampine54(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0 0 0 0 0
               1//6 1//6 0 0 0 0 0 0
               2//9 2//27 4//27 0 0 0 0 0
               3//7 183//1372 -162//343 1053//1372 0 0 0 0
               2//3 68//297 -4//11 42//143 1960//3861 0 0 0
               3//4 597//22528 81//352 63099//585728 58653//366080 4617//20480 0 0
               1 174197//959244 -30942//79937 8152137//19744439 666106//1039181 -29421//29068 482048//414219 0
               1 587//8064 0 4440339//15491840 24353//124800 387//44800 2152//5985 7267//94080
               1 6059//80640 0 8559189//30983680 26411//124800 -927//89600 443//1197 7267//94080]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Bogacki_Shampine_54, butcher)
end

function Tsitouras54(; precision::Type{<:AbstractFloat} = Float64)
    # TEMP until can fix get more digits
    if precision == BigFloat
        @warn "Tsitouras54 can't use BigFloat right now (default to Float64)"
        precision = Float64
    end

    butcher = [0 0 0 0 0 0 0 0
               0.161 0.161 0 0 0 0 0 0
               0.327 -8.4806554923569887e-3 3.3548065549235700e-1 0 0 0 0 0
               0.9 2.8971530571054935 -6.3594484899750752 4.3622954328695815 0 0 0 0
               9.8002554090450966e-1 5.3258648284392569 -1.1748883564062828e1 7.4955393428898365 -9.2495066361755252e-2 0 0 0
               1 5.8614554429464203 -1.2920969317847110e1 8.1593678985761589 -7.1584973281400996e-2 -2.8269050394068383e-2 0 0
               1 9.6460766818065230e-2 0.01 4.7988965041449960e-1 1.3790085741037419 -3.2900695154360808 2.3247105240997739 0
               1 9.6460766818065230e-2 0.01 4.7988965041449960e-1 1.3790085741037419 -3.2900695154360808 2.3247105240997739 0
               1 9.4680755765839453e-2 9.1835655403432540e-3 4.8777052842476160e-1 1.2342975669304790 -2.7077123499835256 1.8666284181705870 1//66]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Tsitouras_54, butcher)
end

function Verner56(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0 0 0 0 0 0
               1//18 1//18 0 0 0 0 0 0 0
               1//6 -1//12 1//4 0 0 0 0 0 0
               2//9 -2//81 4//27 8//81 0 0 0 0 0
               2//3 40//33 -4//11 -56//11 54//11 0 0 0 0
               1 -369//73 72//73 5380//219 -12285//584 2695//1752 0 0 0
               8//9 -8716//891 656//297 39520//891 -416//11 52//27 0 0 0
               1 3015//256 -9//4 -4219//78 5985//128 -539//384 0 693//3328 0
               1 3//80 0 4//25 243//1120 77//160 73//700 0 0
               1 57//640 0 -16//65 1377//2240 121//320 0 891//8320 2//35]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Verner_56, butcher)
end

function Verner65(; precision::Type{<:AbstractFloat} = Float64)
    butcher = [0 0 0 0 0 0 0 0 0
               1//6 1//6 0 0 0 0 0 0 0
               4//15 4//75 16//75 0 0 0 0 0 0
               2//3 5//6 -8//3 5//2 0 0 0 0 0
               5//6 -165//64 55//6 -425//64 85//96 0 0 0 0
               1 12//5 -8 4015//612 -11//36 88//255 0 0 0
               1//15 -8263//15000 124//75 -643//680 -81//250 2484//10625 0 0 0
               1 3501//1720 -300//43 297275//52632 -319//2322 24068//84065 0 3850//26703 0
               1 3//40 0 875//2244 23//72 264//1955 0 125//11592 43//616
               1 13//160 0 2375//5984 5//16 12//85 3//44 0 0]
    butcher = butcher .|> precision 

    RungeKutta(; name = :Verner_65, butcher)
end
