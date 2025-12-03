
abstract type JacobianVectorMethod end

struct NaiveJacobianVector <: JacobianVectorMethod
    evaluations::MVector{1,Int64}
    time_subroutine::Bool
    runtime::MVector{1,Float64}
end

function NaiveJacobianVector()
    evaluations = MVector{1,Int64}(0)
    time_subroutine = false
    runtime = MVector{1,Float64}(0.0)

    return NaiveJacobianVector(evaluations, time_subroutine, runtime)
end

# TODO: check if any T-type bugs
struct FiniteJacobianVector{T} <: JacobianVectorMethod where T <: AbstractFloat
    cache_1::Vector{T}
    cache_2::Vector{T}
    evaluations::MVector{1,Int64}
    time_subroutine::Bool
    runtime::MVector{1,Float64}
end

function FiniteJacobianVector()
    cache_1 = [0.0]
    cache_2 = [0.0]
    evaluations = MVector{1,Int64}(0)
    time_subroutine = false
    runtime = MVector{1,Float64}(0.0)

    return FiniteJacobianVector(cache_1, cache_2, evaluations, time_subroutine, runtime)
end

struct ForwardJacobianVector{D1, D2} <: JacobianVectorMethod where {D1 <: Dual{DeivVecTag},
                                                                    D2 <: Dual{DeivVecTag}}
    cache_1::Vector{D1}
    cache_2::Vector{D2}
    evaluations::MVector{1,Int64}
    time_subroutine::Bool
    runtime::MVector{1,Float64}
end

function ForwardJacobianVector()
    cache_1 = Dual{DeivVecTag}.([0.0], [0.0])
    cache_2 = Dual{DeivVecTag}.([0.0], [0.0])
    evaluations = MVector{1,Int64}(0)
    time_subroutine = false
    runtime = MVector{1,Float64}(0.0)

    return ForwardJacobianVector(cache_1, cache_2, evaluations, time_subroutine, runtime)
end

function reconstruct_jacobian_vector(::NaiveJacobianVector, f::Vector{T},
                                     time_subroutine::Bool) where T <: AbstractFloat
    evaluations = MVector{1,Int64}(0)
    runtime = MVector{1,Float64}(0.0)

    return NaiveJacobianVector(evaluations, time_subroutine, runtime)
end

function reconstruct_jacobian_vector(::FiniteJacobianVector, f::Vector{T},
                                     time_subroutine::Bool) where T <: AbstractFloat
    cache_1 = similar(f)
    cache_2 = similar(f)
    evaluations = MVector{1,Int64}(0)
    runtime = MVector{1,Float64}(0.0)

    return FiniteJacobianVector(cache_1, cache_2, evaluations, time_subroutine, runtime)
end

function reconstruct_jacobian_vector(::ForwardJacobianVector, f::Vector{T},
                                     time_subroutine::Bool) where T <: AbstractFloat
    cache_1 = Dual{DeivVecTag}.(f, f)
    cache_2 = Dual{DeivVecTag}.(f, f)
    evaluations = MVector{1,Int64}(0)
    runtime = MVector{1,Float64}(0.0)

    return ForwardJacobianVector(cache_1, cache_2, evaluations, time_subroutine, runtime)
end

function evaluate_jacobian_vector!(jacobian_vector::NaiveJacobianVector,
                                   JS::SubArray{T}, ode_wrap!::ODEWrapperState,
                                   state_jacobian::JacobianMethod,
                                   J::Union{Matrix{T}, SparseMatrixCSC{T,Int64}},
                                   S::Matrix{T}, y::Vector{T},
                                   f::Vector{T}) where T <: AbstractFloat

    evaluations = jacobian_vector.evaluations
    time_subroutine = jacobian_vector.time_subroutine
    runtime = jacobian_vector.runtime

    runtime_Jy = state_jacobian.runtime[1]

    if time_subroutine && evaluations[1] % SAMPLE_INTERVAL == 0
        stats = @timed begin
            # note: still allocate J for sensitivity methods that don't use it directly
            evaluate_jacobian!(state_jacobian, J, ode_wrap!, y, f)

            # only jacobian-vector routine should time this jacobian evaluation
            state_jacobian.runtime[1] = runtime_Jy
            state_jacobian.evaluations[1] -= 1

            # runtime isn't that bad if J is sparse
            mul!(JS, J, S)
        end
        runtime[1] += SAMPLE_INTERVAL*stats.time
    else
        evaluate_jacobian!(state_jacobian, J, ode_wrap!, y, f)
        state_jacobian.runtime[1] = runtime_Jy
        state_jacobian.evaluations[1] -= 1
        mul!(JS, J, S)
    end
    evaluations[1] += 1

    return nothing
end

function evaluate_jacobian_vector!(jacobian_vector::FiniteJacobianVector,
                                   JS::SubArray{T}, ode_wrap!::ODEWrapperState,
                                   state_jacobian::JacobianMethod,
                                   J::Union{Matrix{T}, SparseMatrixCSC{T,Int64}},
                                   S::Matrix{T}, y::Vector{T},
                                   f::Vector{T}) where T <: AbstractFloat

    cache_1 = jacobian_vector.cache_1
    cache_2 = jacobian_vector.cache_2
    evaluations = jacobian_vector.evaluations
    time_subroutine = jacobian_vector.time_subroutine
    runtime = jacobian_vector.runtime

    p = ode_wrap!.p

    if time_subroutine && evaluations[1] % SAMPLE_INTERVAL == 0
        stats = @timed begin
            for j in eachindex(p)
                Jv = view(JS, :, j)
                v = view(S, :, j)
                num_jacvec_tmp!(Jv, ode_wrap!, y, v, f, cache_1, cache_2)
            end
        end
        runtime[1] += SAMPLE_INTERVAL*stats.time
    else
        for j in eachindex(p)
            Jv = view(JS, :, j)
            v = view(S, :, j)
            num_jacvec_tmp!(Jv, ode_wrap!, y, v, f, cache_1, cache_2)
        end
    end
    evaluations[1] += 1

    return nothing
end

function evaluate_jacobian_vector!(jacobian_vector::ForwardJacobianVector,
                                   JS::SubArray{T}, ode_wrap!::ODEWrapperState,
                                   state_jacobian::JacobianMethod,
                                   J::Union{Matrix{T}, SparseMatrixCSC{T,Int64}},
                                   S::Matrix{T}, y::Vector{T},
                                   args...) where T <: AbstractFloat

    cache_1 = jacobian_vector.cache_1
    cache_2 = jacobian_vector.cache_2
    evaluations = jacobian_vector.evaluations
    time_subroutine = jacobian_vector.time_subroutine
    runtime = jacobian_vector.runtime

    p = ode_wrap!.p

    if time_subroutine && evaluations[1] % SAMPLE_INTERVAL == 0
        stats = @timed begin
            for j in eachindex(p)
                Jv = view(JS, :, j)
                v = view(S, :, j)
                auto_jacvec!(Jv, ode_wrap!, y, v, cache_1, cache_2)
            end
        end
        runtime[1] += SAMPLE_INTERVAL*stats.time
    else
        for j in eachindex(p)
            Jv = view(JS, :, j)
            v = view(S, :, j)
            auto_jacvec!(Jv, ode_wrap!, y, v, cache_1, cache_2)
        end
    end
    evaluations[1] += 1

    return nothing
end

# TODO: override method from SparseDiffTools
# note: modified version of num_jacvec! from SparseDiffTools
function num_jacvec_tmp!(Jv, ode_wrap!, y, v, f, cache_1, cache_2;
                         epsilon = sqrt(eps(1.0)), alpha = epsilon)
    v_norm = sqrt(sum(abs2, v))

    if iszero(v_norm)
        @.. Jv = 0.0
    else
        lambda = max(alpha, epsilon*abs(dot(y, v))/v_norm)
        @.. cache_2 = y + lambda*v
        ode_wrap!(cache_1, cache_2)
        @.. Jv = (cache_1 - f) / lambda
    end
    return nothing
end
