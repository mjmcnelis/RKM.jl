
# TODO: get more digits and try BigFloat
function Feagin108(; precision::Type{<:AbstractFloat} = Float64)
    if precision == BigFloat || precision == Double64
        @warn "Feagin108 can't use $precision right now (default to Float64)"
        precision = Float64
    end
    butcher = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
               1//10 1//10 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
               53935784080298178753//100000000000000000000 -22879414034382286013//25000000000000000000 29090688043565464561//20000000000000000000 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
               8090367612044726813//10000000000000000000 2528239878763977129//12500000000000000000 0 60677757090335451097//100000000000000000000 0 0 0 0 0 0 0 0 0 0 0 0 0 0
               3090367612044726813//10000000000000000000 3680494294172871503//20000000000000000000 0 19796683122719236907//100000000000000000000 -1823869618284081573//25000000000000000000 0 0 0 0 0 0 0 0 0 0 0 0 0
               3924296760879181073//4000000000000000000 8790073402066813373//100000000000000000000 0 0 10261492563006516133//25000000000000000000 1206784384197166223//2500000000000000000 0 0 0 0 0 0 0 0 0 0 0 0
               83333333333333333333//100000000000000000000 4298502524512301511//50000000000000000000 0 0 6617719260814443679//20000000000000000000 12241573932736254821//25000000000000000000 -914820468838563421//12500000000000000000 0 0 0 0 0 0 0 0 0 0 0
               35401736585680237633//100000000000000000000 6046522456266686033//50000000000000000000 0 0 0 26012467575829562281//100000000000000000000 1627013107745456651//50000000000000000000 -186181316192925313//3125000000000000000 0 0 0 0 0 0 0 0 0 0
               88252766196473234643//100000000000000000000 11085437958039148351//100000000000000000000 0 0 0 0 -1514403720637513969//25000000000000000000 3217637056017783901//10000000000000000000 25524286280403151579//50000000000000000000 0 0 0 0 0 0 0 0 0
               2008174244501007963//3125000000000000000 11205441475287900483//100000000000000000000 0 0 0 0 -14494277590286591567//100000000000000000000 -33326971909625670659//100000000000000000000 9985384591137601227//20000000000000000000 6368807611621076303//12500000000000000000 0 0 0 0 0 0 0 0
               1116825755498992037//3125000000000000000 5698839198209299307//50000000000000000000 0 0 0 0 -3844066821016784693//50000000000000000000 23952736032439064911//100000000000000000000 7955493247361892781//20000000000000000000 268897392184018639//25000000000000000000 -6555382483280377483//20000000000000000000 0 0 0 0 0 0 0
               11747233803526765357//100000000000000000000 498946580175122529//6250000000000000000 0 0 0 0 -1040659373601206153//20000000000000000000 -5769541461685488817//100000000000000000000 9739095785605208249//50000000000000000000 14538492318832506973//100000000000000000000 -122334798492448559//1562500000000000000 -5725164968054945609//50000000000000000000 0 0 0 0 0 0
               83333333333333333333//100000000000000000000 24627890254121432003//25000000000000000000 0 0 6617719260814443679//20000000000000000000 12241573932736254821//25000000000000000000 -68948243287421783791//50000000000000000000 -86116419502763566667//100000000000000000000 578428813637537220023//100000000000000000000 32880776198510356689//10000000000000000000 -238633905093136384013//100000000000000000000 -65095868496728783731//20000000000000000000 -108171770843211491177//50000000000000000000 0 0 0 0 0
               3090367612044726813//10000000000000000000 17901605915432657821//20000000000000000000 0 19796683122719236907//100000000000000000000 -1823869618284081573//25000000000000000000 0 -42561811983100380987//50000000000000000000 9958002807963332543//25000000000000000000 363937263181035606029//100000000000000000000 30964575407966064473//20000000000000000000 -106110857352026858013//50000000000000000000 -158350398545326172713//100000000000000000000 -85780804142968132461//50000000000000000000 -2440364057501274521//100000000000000000000 0 0 0 0
               53935784080298178753//100000000000000000000 -22879414034382286013//25000000000000000000 29090688043565464561//20000000000000000000 0 0 -38866682182248411677//50000000000000000000 0 -1138619577693970087//12500000000000000000 0 0 0 0 0 1138619577693970087//12500000000000000000 38866682182248411677//50000000000000000000 0 0 0
               1//10 1//10 0 -15717866579977116337//100000000000000000000 0 0 0 0 0 0 0 0 0 0 0 15717866579977116337//100000000000000000000 0 0
               1 18178130070009528389//100000000000000000000 27//40 17137907992359491997//50000000000000000000 0 25911121454832274451//100000000000000000000 -7165579334359041781//20000000000000000000 -20918979188176661219//20000000000000000000 93032784541562698329//100000000000000000000 88975479715854051223//50000000000000000000 1//10 -28254756953904408161//100000000000000000000 -15932735011997254917//100000000000000000000 -7275794732350075543//50000000000000000000 -25911121454832274451//100000000000000000000 -17137907992359491997//50000000000000000000 -27//40 0
               1 1//30 1//40 1//30 0 1//20 0 1//25 0 2365468476861543627//12500000000000000000 27742918851774317651//100000000000000000000 27742918851774317651//100000000000000000000 2365468476861543627//12500000000000000000 -1//25 -1//20 -1//30 -1//40 1//30
               1 1//30 1//36 1//30 0 1//20 0 1//25 0 2365468476861543627//12500000000000000000 27742918851774317651//100000000000000000000 27742918851774317651//100000000000000000000 2365468476861543627//12500000000000000000 -1//25 -1//20 -1//30 -1//36 1//30]
    butcher = butcher .|> precision

    RungeKutta(; name = :Feagin_10_8, butcher)
end
