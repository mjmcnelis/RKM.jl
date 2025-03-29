
abstract type JacobianVectorMethod end

# TODO: consider renaming this to AutoJacVec
@kwdef struct ForwardJacobianVector{D1, D2} <: JacobianVectorMethod where {
                                                     D1 <: Dual{DeivVecTag},
                                                     D2 <: Dual{DeivVecTag}}
    cache_1::Vector{D1} = Dual{DeivVecTag}.([0.0], [0.0])
    cache_2::Vector{D2} = Dual{DeivVecTag}.([0.0], [0.0])
end

# TODO: check if any T-type bugs
@kwdef struct FiniteJacobianVector{T} <: JacobianVectorMethod where T <: AbstractFloat
    cache_1::Vector{T} = [0.0]
    cache_2::Vector{T} = [0.0]
end