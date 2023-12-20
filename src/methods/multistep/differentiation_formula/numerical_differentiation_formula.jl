"""
    NumericalDifferentiationFormula(; order::Int64,
        precision::Type{T} = Float64) where T <: AbstractFloat

Numerical differentiation formula (NDF) implicit multistep method.

Note: `order` ranges from 1-4, `table_pred` contains the predictor coefficients
"""
function NumericalDifferentiationFormula(; order::Int64,
            precision::Type{T} = Float64) where T <: AbstractFloat

    @assert 1 <= order <= 4 "order = $order is not valid"

    κ = [-37//200, -1//9, -823//10000, -83//2000]
    γ = [1, 3//2, 11//6, 25//12]
    κγ = κ .* γ

    κγ1 = κγ[1]
    κγ2 = κγ[2]
    κγ3 = κγ[3]
    κγ4 = κγ[4]

    # TODO: as it stands, DifferentiationFormula(; ...) won't pick up last nonzero coefficient
    table = [1/(1 - κγ1) (1 - 2κγ1)/(1 - κγ1) κγ1/(1 - κγ1) 0 0 0
             2//3/(1 - κγ2) (4//3 - 3κγ2)/(1 - κγ2) -(1//3 - 3κγ2)/(1 - κγ2) -κγ2/(1 - κγ2) 0 0
             6//11/(1 - κγ3) (18//11 - 4κγ3)/(1 - κγ3) -(9//11 - 6κγ3)/(1 - κγ3) (2//11 - 4κγ3)/(1 - κγ3) κγ3/(1 - κγ3) 0
             12//25/(1 - κγ4) (48//25 - 5κγ4)/(1 - κγ4) -(36//25 - 10κγ4)/(1 - κγ4) (16//25 - 10κγ4)/(1 - κγ4) -(3//25 - 5κγ4)/(1 - κγ4) -κγ4/(1 - κγ4)
            ]
    table = table .|> precision

    table_pred = [0 1 -1 0 0 0
                  0 2 -3 1 0 0
                  0  3 -6 4 -1 0
                  0 4 -10 10 -5 1
                 ]
    table_pred = table_pred .|> precision

    name = Symbol("Numerical_Differentiation_Formula_$order")

    return DifferentiationFormula(; name, order, table, table_pred)
end
