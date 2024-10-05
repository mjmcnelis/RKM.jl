
abstract type SensitivityMethod end

struct NoSensitivity <: SensitivityMethod end

# is there a better name?
@kwdef struct DecoupledDirect{JM} <: SensitivityMethod
    jacobian_method::JM = FiniteJacobian()
end