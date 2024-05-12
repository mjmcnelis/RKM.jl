
# TODO: break up properties in directory

abstract type Iteration end
abstract type Implicit <: Iteration end

struct Explicit <: Iteration end

# implicit Runge-Kutta methods
struct DiagonalImplicit <: Implicit end
struct FullImplicit <: Implicit end

# implicit multistep methods
struct SingleImplicit <: Implicit end

# TODO: determine if butcher table fixed (square matrix) or embedded (not square matrix)

function method_is_fsal(butcher::SMatrix{N, M, T, NM}) where {N, M, T <: AbstractFloat, NM}
    # TODO: BackwardEuler1 is incorrectly labeled as fsal = true
    #       need to check if explicit stage[1] = true
    ncol = size(butcher, 2)

    # TODO: if using implicit routine TRBDF2, then also need to check
    #       whether first stage is explicit in order to use f = f_tmp
    #       should probably test it out on TRBDF2 (fixed time step and embedded case)

    # remove any embedded rows
    B_square = butcher[1:ncol, :]

    # check if last, second-last rows are identical
    if B_square[end, :] == B_square[end-1, :]
        fsal = true
    else
        fsal = false
    end
    return fsal
end

# note: needed @inline macro to make SVector return type-stable
@inline function order_prop(name::Symbol, precision::Type{T},
                            p::Int) where T <: AbstractFloat
    order = filter.(isdigit, split(string(name), "_"))
    order = parse.(precision, filter(x -> x != "", order))
    return SVector{p, T}(order)
end

function explicit_stage_prop(butcher::SMatrix{N, M, T, NM}) where {N, M,
                                                                   T <: AbstractFloat,
                                                                   NM}
    ncol = size(butcher, 2)
    A = butcher[1:(ncol-1), 2:end]
    stages = ncol-1
    explicit_stage = Vector{Bool}(undef, stages)
    for i in 1:stages
        explicit_stage[i] = all(x -> x == 0.0, A[i,i:end])
    end
    return explicit_stage |> SVector{stages}
end