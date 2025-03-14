"""
    Fehlberg45(precision::Type{T} = Float64) where T <: AbstractFloat

Fehlberg's fourth(fifth)-order method.

https://ntrs.nasa.gov/citations/19690021375
"""
function Fehlberg45(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Fehlberg_4_5
    butcher = SMatrix{7, 8, precision, 56}(
        0, 0, 0, 0, 0, 0, 0,
        1//4, 1//4, 0, 0, 0, 0, 0,
        3//8, 3//32, 9//32, 0, 0, 0, 0,
        12//13, 1932//2197, -7200//2197, 7296//2197, 0, 0, 0,
        1, 439//216, -8, 3680//513, -845//4104, 0, 0,
        1//2, -8//27, 2, -3544//2565, 1859//4104, -11//40, 0,
        1, 25//216, 0, 1408//2565, 2197//4104, -1//5, 0,
        1, 16//135, 0, 6656//12825, 28561//56430, -9//50, 2//55
    ) |> transpose
    iteration = Explicit()
    reconstructor = Fehlberg45

    return RungeKutta(name, butcher, iteration, reconstructor)
end

"""
    CashKarp54(precision::Type{T} = Float64) where T <: AbstractFloat

Cash and Karp's fifth(fourth)-order method.

http://www.elegio.it/mc2/rk/doc/p201-cash-karp.pdf
"""
function CashKarp54(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Cash_Karp_5_4
    butcher = SMatrix{7, 8, precision, 56}(
        0, 0, 0, 0, 0, 0, 0,
        1//5, 1//5, 0, 0, 0, 0, 0,
        3//10, 3//40, 9//40, 0, 0, 0, 0,
        3//5, 3//10, -9//10, 6//5, 0, 0, 0,
        1, -11//54, 5//2, -70//27, 35//27, 0, 0,
        7//8, 1631//55296, 175//512, 575//13824, 44275//110592, 253//4096, 0,
        1, 37//378, 0, 250//621, 125//594, 0, 512//1771,
        1, 2825//27648, 0, 18575//48384, 13525//55296, 277//14336, 1//4
    ) |> transpose
    iteration = Explicit()
    reconstructor = CashKarp54

    return RungeKutta(name, butcher, iteration, reconstructor)
end

"""
    DormandPrince54(precision::Type{T} = Float64) where T <: AbstractFloat

Dormand and Prince's fifth(fourth)-order method.

https://www.sciencedirect.com/science/article/pii/0771050X80900133?via%3Dihub
"""
function DormandPrince54(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Dormand_Prince_5_4
    butcher = SMatrix{8, 9, precision, 72}(
        0, 0, 0, 0, 0, 0, 0, 0,
        1//5, 1//5, 0, 0, 0, 0, 0, 0,
        3//10, 3//40, 9//40, 0, 0, 0, 0, 0,
        4//5, 44//45, -56//15, 32//9, 0, 0, 0, 0,
        8//9, 19372//6561, -25360//2187, 64448//6561, -212//729, 0, 0, 0,
        1, 9017//3168, -355//33, 46732//5247, 49//176, -5103//18656, 0, 0,
        1, 35//384, 0, 500//1113, 125//192, -2187//6784, 11//84, 0,
        1, 35//384, 0, 500//1113, 125//192, -2187//6784, 11//84, 0,
        1, 5179//57600, 0, 7571//16695, 393//640, -92097//339200, 187//2100, 1//40
    ) |> transpose
    iteration = Explicit()
    reconstructor = DormandPrince54

    return RungeKutta(name, butcher, iteration, reconstructor)
end

"""
    BogackiShampine54(precision::Type{T} = Float64) where T <: AbstractFloat

Bogacki and Shampine's fifth(fourth)-order method.
"""
function BogackiShampine54(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Bogacki_Shampine_5_4
    butcher = SMatrix{8, 9, precision, 72}(
        0, 0, 0, 0, 0, 0, 0, 0,
        1//6, 1//6, 0, 0, 0, 0, 0, 0,
        2//9, 2//27, 4//27, 0, 0, 0, 0, 0,
        3//7, 183//1372, -162//343, 1053//1372, 0, 0, 0, 0,
        2//3, 68//297, -4//11, 42//143, 1960//3861, 0, 0, 0,
        3//4, 597//22528, 81//352, 63099//585728, 58653//366080, 4617//20480, 0, 0,
        1, 174197//959244, -30942//79937, 8152137//19744439, 666106//1039181, -29421//29068, 482048//414219, 0,
        1, 587//8064, 0, 4440339//15491840, 24353//124800, 387//44800, 2152//5985, 7267//94080,
        1, 6059//80640, 0, 8559189//30983680, 26411//124800, -927//89600, 443//1197, 7267//94080
    ) |> transpose
    iteration = Explicit()
    reconstructor = BogackiShampine54

    return RungeKutta(name, butcher, iteration, reconstructor)
end

"""
    Tsitouras54(precision::Type{T} = Float64) where T <: AbstractFloat

Tsitouras' fifth(fourth)-order method.

https://www.sciencedirect.com/science/article/pii/S0898122111004706
"""
function Tsitouras54(precision::Type{T} = Float64) where T <: AbstractFloat
    # TODO: not sure why other tilde pair (commented row) sums to (0.97, paper)
    #       or (0, OrdinaryDiffEq) instead of 1
    name = :Tsitouras_5_4
    butcher = SMatrix{8, 9, precision, 72}(
        0, 0, 0, 0, 0, 0, 0, 0,
        161//1000, 161//1000, 0, 0, 0, 0, 0, 0,
        327//1000, -big"8.480655492356988544426874250230774675121177393430391537369234245294192976164099e-3", big"3.354806554923569885444268742502307746751211773934303915373692342452941929761659e-1", 0, 0, 0, 0, 0,
        9//10, big"2.897153057105493432130432594192938764924887287701866490314866693455023795137512", -big"6.359448489975074843148159912383825625952700647415626703305928850207288721235239", big"4.362295432869581411017727318190886861027813359713760212991062156752264926097707", 0, 0, 0, 0,
        big"9.800255409045096857298102862870245954942137979563024768854764293221195950761096e-1", big"5.325864828439256604428877920840511317836476253097040101202360397727981648835585", -big"11.74888356406282787774717033978577296188744178259862899288666928009020615663589", big"7.495539342889836208304604784564358155658679161518186721010132816213648793440556", -big"9.24950663617552492565020793320719161134998340602953524403475045293046905641141e-2", 0, 0, 0,
        1, big"5.861455442946420028659251486982647890394337666164814434818157239052507339770685", -big"1.292096931784710929170611868178335939541780751955743459166312250439928519268342e1", big"8.159367898576158643180400794539253485181918321135053305748355423955009222648607", -big"7.158497328140099722453054252582973869127213147363544882721139659546372402303727e-2", -big"2.826905039406838290900305721271224146717633626879770007617876201276764571291576e-2", 0, 0,
        1, big"9.646076681806522951816731316512876333711995238157997181903319145764851595234052e-2", 1//100, big"4.798896504144995747752495322905965199130404621990332488332634944254542060153088e-1", big"1.379008574103741893192274821856872770756462643091360525934940067397245698027568", -big"3.290069515436080679901047585711363850115683290894936158531296799594813811049913", big"2.324710524099773982415355918398765796109060233222962411944060046314465391054717", 0,
        1, big"9.646076681806522951816731316512876333711995238157997181903319145764851595234052e-2", 1//100, big"4.798896504144995747752495322905965199130404621990332488332634944254542060153088e-1", big"1.379008574103741893192274821856872770756462643091360525934940067397245698027568", -big"3.290069515436080679901047585711363850115683290894936158531296799594813811049913", big"2.324710524099773982415355918398765796109060233222962411944060046314465391054717", 0,
        #1, big"1.780011052225771443378550607539534775944678804333659557637450799792588061629796e-3", big"8.164344596567469032236360633546862401862537590159047610940604670770447527463931e-4", -big"7.880878010261996010314727672526304238628733777103128603258129604952959142646516e-3", big"1.44711007173262907537165147972635116720922712343167677619514233896760819649515e-1", -big"5.823571654525552250199376106520421794260781239567387797673045438803694038950012e-1", big"4.580821059291869466616365188325542974428047279788398179474684434732070620889539e-1", 1//66,
        1, big"9.468075576583945807478876255758922856117527357724631226139574065785592789071072e-2", big"9.183565540343253096776363936645313759813746240984095238905939532922955247253601e-3", big"4.877705284247615707855642599631228241516691959761363774365216240304071651579553e-1", big"1.234297566930478985655109673884237654035539930748192848315425833500484878378053", -big"2.707712349983525454881109975059321670689605166938197378763992255714444407154911", big"1.866628418170587035753719399566211498666255505244122593996591602841258328965763", 1//66
    ) |> transpose

    # polynomial coefficients for continuous output
    ω = SMatrix{4, 7, Float64, 28}(
        0.9999999999999999, -2.763706197274826, 2.9132554618219126, -1.0530884977290216,
        0.0, 0.13169999999999998, -0.2234, 0.1017,
        0.0, 3.930296236894751, -5.941033872131505, 2.490627285651253,
        0.0, -12.411077166933676, 30.338188630282318, -16.548102889244902,
        0.0, 37.50931341651104, -88.1789048947664, 47.37952196281928,
        0.0, -27.896526289197286, 65.09189467479368, -34.87065786149661,
        0.0, 1.5, -4.0, 2.5
    ) |> transpose

    iteration = Explicit()
    reconstructor = Tsitouras54

    # return RungeKutta(name, butcher, ω, iteration, reconstructor)
    return RungeKutta(name, butcher, iteration, reconstructor; ω)
end

"""
    Verner56(precision::Type{T} = Float64) where T <: AbstractFloat

Verner's fifth(sixth)-order method (1978).

https://www.jstor.org/stable/2156853
"""
function Verner56(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Verner_5_6
    butcher = SMatrix{9, 10, precision, 90}(
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        1//18, 1//18, 0, 0, 0, 0, 0, 0, 0,
        1//6, -1//12, 1//4, 0, 0, 0, 0, 0, 0,
        2//9, -2//81, 4//27, 8//81, 0, 0, 0, 0, 0,
        2//3, 40//33, -4//11, -56//11, 54//11, 0, 0, 0, 0,
        1, -369//73, 72//73, 5380//219, -12285//584, 2695//1752, 0, 0, 0,
        8//9, -8716//891, 656//297, 39520//891, -416//11, 52//27, 0, 0, 0,
        1, 3015//256, -9//4, -4219//78, 5985//128, -539//384, 0, 693//3328, 0,
        1, 3//80, 0, 4//25, 243//1120, 77//160, 73//700, 0, 0,
        1, 57//640, 0, -16//65, 1377//2240, 121//320, 0, 891//8320, 2//35
    ) |> transpose
    iteration = Explicit()
    reconstructor = Verner56

    return RungeKutta(name, butcher, iteration, reconstructor)
end

"""
    Verner65(precision::Type{T} = Float64) where T <: AbstractFloat

Verner's sixth(fifth)-order method.

https://link.springer.com/book/10.1007/978-3-540-78862-1
"""
function Verner65(precision::Type{T} = Float64) where T <: AbstractFloat
    name = :Verner_6_5
    butcher = SMatrix{9, 10, precision, 90}(
        0, 0, 0, 0, 0, 0, 0, 0, 0,
        1//6, 1//6, 0, 0, 0, 0, 0, 0, 0,
        4//15, 4//75, 16//75, 0, 0, 0, 0, 0, 0,
        2//3, 5//6, -8//3, 5//2, 0, 0, 0, 0, 0,
        5//6, -165//64, 55//6, -425//64, 85//96, 0, 0, 0, 0,
        1, 12//5, -8, 4015//612, -11//36, 88//255, 0, 0, 0,
        1//15, -8263//15000, 124//75, -643//680, -81//250, 2484//10625, 0, 0, 0,
        1, 3501//1720, -300//43, 297275//52632, -319//2322, 24068//84065, 0, 3850//26703, 0,
        1, 3//40, 0, 875//2244, 23//72, 264//1955, 0, 125//11592, 43//616,
        1, 13//160, 0, 2375//5984, 5//16, 12//85, 3//44, 0, 0
    ) |> transpose
    iteration = Explicit()
    reconstructor = Verner65

    return RungeKutta(name, butcher, iteration, reconstructor)
end
