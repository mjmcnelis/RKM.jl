
abstract type JacobianVectorMethod end

struct NaiveJacobianVector <: JacobianVectorMethod end

# TODO: check if any T-type bugs
@kwdef struct FiniteJacobianVector{T} <: JacobianVectorMethod where T <: AbstractFloat
    cache_1::Vector{T} = [0.0]
    cache_2::Vector{T} = [0.0]
end

# TODO: consider renaming this to AutoJacVec
@kwdef struct ForwardJacobianVector{D1, D2} <: JacobianVectorMethod where {
                                                     D1 <: Dual{DeivVecTag},
                                                     D2 <: Dual{DeivVecTag}}
    cache_1::Vector{D1} = Dual{DeivVecTag}.([0.0], [0.0])
    cache_2::Vector{D2} = Dual{DeivVecTag}.([0.0], [0.0])
end

function reconstruct_jacobian_vector(::NaiveJacobianVector, args...)
    return NaiveJacobianVector()
end

function reconstruct_jacobian_vector(::FiniteJacobianVector,
                                     f::Vector{T}) where T <: AbstractFloat
    cache_1 = similar(f)
    cache_2 = similar(f)
    return FiniteJacobianVector(; cache_1, cache_2)
end

function reconstruct_jacobian_vector(::ForwardJacobianVector,
                                     f::Vector{T}) where T <: AbstractFloat
    cache_1 = Dual{DeivVecTag}.(f, f)
    cache_2 = Dual{DeivVecTag}.(f, f)
    return ForwardJacobianVector(; cache_1, cache_2)
end

function evaluate_jacobian_sensitivity!(jacobian_vector::NaiveJacobianVector,
                                        JS::SubArray{T}, ode_wrap!::ODEWrapperState,
                                        state_jacobian::JacobianMethod,
                                        J::Union{Matrix{T}, SparseMatrixCSC{T,Int64}},
                                        S::Matrix{T}, y::Vector{T},
                                        f::Vector{T}) where T <: AbstractFloat
    # note: still allocate J for sensitivity methods that don't use it directly
    evaluate_jacobian!(state_jacobian, J, ode_wrap!, y, f)
    # runtime isn't that bad if J is sparse
    mul!(JS, J, S)
    return nothing
end

function evaluate_jacobian_sensitivity!(jacobian_vector::FiniteJacobianVector,
                                        JS::SubArray{T}, ode_wrap!::ODEWrapperState,
                                        state_jacobian::JacobianMethod,
                                        J::Union{Matrix{T}, SparseMatrixCSC{T,Int64}},
                                        S::Matrix{T}, y::Vector{T},
                                        f::Vector{T}) where T <: AbstractFloat

    cache_1 = jacobian_vector.cache_1
    cache_2 = jacobian_vector.cache_2

    p = ode_wrap!.p

    # TODO: type-dispatch Jv subroutine
    for j in eachindex(p)
        Jv = view(JS, :, j)
        v = view(S, :, j)
        num_jacvec_tmp!(Jv, ode_wrap!, y, v, f, cache_1, cache_2)
    end
    return nothing
end

function evaluate_jacobian_sensitivity!(jacobian_vector::ForwardJacobianVector,
                                        JS::SubArray{T}, ode_wrap!::ODEWrapperState,
                                        state_jacobian::JacobianMethod,
                                        J::Union{Matrix{T}, SparseMatrixCSC{T,Int64}},
                                        S::Matrix{T}, y::Vector{T},
                                        args...) where T <: AbstractFloat

    cache_1 = jacobian_vector.cache_1
    cache_2 = jacobian_vector.cache_2

    p = ode_wrap!.p

    for j in eachindex(p)
        Jv = view(JS, :, j)
        v = view(S, :, j)
        auto_jacvec!(Jv, ode_wrap!, y, v, cache_1, cache_2)
    end
    return nothing
end

# TODO: override method from SparseDiffTools
# note: modified version of num_jacvec! from SparseDiffTools
function num_jacvec_tmp!(Jv, ode_wrap!, y, v, f, cache_1, cache_2;
                         epsilon = sqrt(eps(1.0)), alpha = epsilon)
    v_norm = sqrt(sum(abs2, v))

    if v_norm == 0.0
        @.. Jv = 0.0
    else
        lambda = max(alpha, epsilon*abs(dot(y, v))/v_norm)
        @.. cache_2 = y + lambda*v
        ode_wrap!(cache_1, cache_2)
        @.. Jv = (cache_1 - f) / lambda
    end
    return nothing
end
