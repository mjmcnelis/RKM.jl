"""
    BackwardDifferentiationFormula(;
        precision::Type{T} = Float64) where T <: AbstractFloat

Backward differentiation formula (BDF) implicit multistep method.

Note: `order` ranges from 1-6, `table_pred` contains the predictor coefficients
"""
function BackwardDifferentiationFormula(; order::Int64,
            precision::Type{T} = Float64) where T <: AbstractFloat

    @assert 1 <= order <= 6 "order = $order is not valid"

    table = [1 1 0 0 0 0 0
             2//3 4//3 -1//3 0 0 0 0
             6//11 18//11 -9//11 2//11 0 0 0
             12//25 48//25 -36//25 16//25 -3//25 0 0
             60//137 300//137 -300//137 200//137 -75//137 12//137 0
             60//147 360//147 -450//147 400//147 -225//147 72//147 -10//147
            ]
    table = table .|> precision

    table_pred = [1 -1 0 0 0 0 0
                  2 -3 1 0 0 0 0
                  3 -6 4 -1 0 0 0
                  4 -10 10 -5 1 0 0
                  5 -15 20 -15 6 -1 0
                  6 -21 35 -35 21 -7 1
                 ]
    table_pred = table_pred .|> precision

    name = Symbol("Backward_Differentiation_Formula_$order")

    return DifferentiationFormula(; name, order, table, table_pred)
end
