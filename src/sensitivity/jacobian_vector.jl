
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

function evaluate_jacobian_sensitivity!(jacobian_vector::FiniteJacobianVector,
                                        JS::SubArray{T}, ode_wrap!::ODEWrapperState,
                                        t::T, S::Matrix{T}, y::Vector{T},
                                        f::Vector{T}) where T <: AbstractFloat
    @unpack cache_1, cache_2 = jacobian_vector
    @unpack p = ode_wrap!
    ode_wrap!.t[1] = t

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
                                        t::T, S::Matrix{T}, y::Vector{T},
                                        args...) where T <: AbstractFloat
    @unpack cache_1, cache_2 = jacobian_vector
    @unpack p = ode_wrap!
    ode_wrap!.t[1] = t

    for j in eachindex(p)
        Jv = view(JS, :, j)
        v = view(S, :, j)
        auto_jacvec!(Jv, ode_wrap!, y, v, cache_1, cache_2)
    end
    return nothing
end

# note: modified version of num_jacvec! from SparseDiffTools
function num_jacvec_tmp!(Jv, ode_wrap!, y, v, f, cache_1, cache_2;
                         epsilon = sqrt(eps(1.0)), alpha = epsilon)
    v_norm = sqrt(sum(abs2, v))

    if v_norm == 0.0
        @.. Jv = 0.0
    else
        λ = max(alpha, epsilon*abs(dot(y, v))/v_norm)
        @.. cache_2 = y + λ*v
        ode_wrap!(cache_1, cache_2)
        @.. Jv = (cache_1 - f) / λ
    end
    return nothing
end
