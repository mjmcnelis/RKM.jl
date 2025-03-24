
abstract type JacobianVectorMethod end

@kwdef struct ForwardJacobianVector{DC1, DC2} <: JacobianVectorMethod where {
                                                     DC1 <: DerivativeConfig,
                                                     DC2 <: DerivativeConfig}
    dcache_1::DC1 = DerivativeConfig(nothing, [0.0], 0.0)
    dcache_2::DC2 = DerivativeConfig(nothing, [0.0], 0.0)
end
