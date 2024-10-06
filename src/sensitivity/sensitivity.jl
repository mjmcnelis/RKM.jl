
abstract type SensitivityMethod end

struct NoSensitivity <: SensitivityMethod end

@kwdef struct DecoupledDirect{JM} <: SensitivityMethod
    jacobian_method::JM = FiniteJacobian()
end