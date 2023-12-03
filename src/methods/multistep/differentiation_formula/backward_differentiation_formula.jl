"""
    BackwardDifferentiationFormula1(;
        precision::Type{T} = Float64) where T <: AbstractFloat

First-order backward differentiation formula or BDF1.

Note: second row in `table` are predictor coefficients
"""
function BackwardDifferentiationFormula1(;
            precision::Type{T} = Float64) where T <: AbstractFloat

    table = [1 1
             1 -1]
    table = table .|> precision
    return DifferentiationFormula(; name = :Backward_Differentiation_Formula_1, table)
end

"""
    BackwardDifferentiationFormula2(;
        precision::Type{T} = Float64) where T <: AbstractFloat

Second-order backward differentiation formula or BDF2.

Note: second row in `table` are predictor coefficients
"""
function BackwardDifferentiationFormula2(;
            precision::Type{T} = Float64) where T <: AbstractFloat

    table = [2//3 4//3 -1//3
             2 -3 1]
    table = table .|> precision
    return DifferentiationFormula(; name = :Backward_Differentiation_Formula_2, table)
end

"""
    BackwardDifferentiationFormula3(;
        precision::Type{T} = Float64) where T <: AbstractFloat

Third-order backward differentiation formula or BDF3.

Note: second row in `table` are predictor coefficients
"""
function BackwardDifferentiationFormula3(;
            precision::Type{T} = Float64) where T <: AbstractFloat

    table = [6//11 18//11 -9//11 2//11
             3 -6 4 -1]
    table = table .|> precision
    return DifferentiationFormula(; name = :Backward_Differentiation_Formula_3, table)
end

"""
    BackwardDifferentiationFormula4(;
        precision::Type{T} = Float64) where T <: AbstractFloat

Fourth-order backward differentiation formula or BDF4.

Note: second row in `table` are predictor coefficients
"""
function BackwardDifferentiationFormula4(;
            precision::Type{T} = Float64) where T <: AbstractFloat

    table = [12//25 48//25 -36//25 16//25 -3//25
             4 -10 10 -5 1]
    table = table .|> precision
    return DifferentiationFormula(; name = :Backward_Differentiation_Formula_4, table)
end

"""
    BackwardDifferentiationFormula5(;
        precision::Type{T} = Float64) where T <: AbstractFloat

Fifth-order backward differentiation formula or BDF5.

Note: second row in `table` are predictor coefficients
"""
function BackwardDifferentiationFormula5(;
            precision::Type{T} = Float64) where T <: AbstractFloat

    table = [60//137 300//137 -300//137 200//137 -75//137 12//137
             5 -15 20 -15 6 -1]
    table = table .|> precision
    return DifferentiationFormula(; name = :Backward_Differentiation_Formula_5, table)
end

"""
    BackwardDifferentiationFormula6(;
        precision::Type{T} = Float64) where T <: AbstractFloat

Sixth-order backward differentiation formula or BDF6.

Note: second row in `table` are predictor coefficients
"""
function BackwardDifferentiationFormula6(;
            precision::Type{T} = Float64) where T <: AbstractFloat

    table = [60//147 360//147 -450//147 400//147 -225//147 72//147 -10//147
             6 -21 35 -35 21 -7 1]
    table = table .|> precision
    return DifferentiationFormula(; name = :Backward_Differentiation_Formula_6, table)
end
