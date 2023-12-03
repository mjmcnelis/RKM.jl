"""
    BackwardDifferentiationFormula1(;
        precision::Type{T} = Float64) where T <: AbstractFloat

First-order backward differentiation formula (BDF1).

Note: second row in `table` are predictor coefficients
"""
function BackwardDifferentiationFormula1(;
            precision::Type{T} = Float64) where T <: AbstractFloat

    table = [1 1
             1 1]
    table = table .|> precision
    return DifferentiationFormula(; name = :Backward_Differentiation_Formula_1, table)
end
