
# TODO: break up properties in directory

abstract type Iteration end
abstract type Implicit <: Iteration end

struct Explicit <: Iteration end

# implicit Runge-Kutta methods
struct DiagonalImplicit <: Implicit end
struct FullImplicit <: Implicit end

# implicit multistep methods
struct SingleImplicit <: Implicit end

abstract type FirstSameAsLast end
struct FSAL <: FirstSameAsLast end
struct NotFSAL <: FirstSameAsLast end

# TODO: determine if butcher table fixed (square matrix) or embedded (not square matrix)

function iteration_prop(butcher::Matrix{T}) where T <: AbstractFloat
    ncol = size(butcher, 2)
    A = butcher[1:(ncol-1), 2:end]

    # check if submatrix A_{ij} is lower triangular
    if tril(A) == A
        # check if diagonal elements of A are all zero
        if all(x -> x == 0, diag(A))
            iteration = Explicit()
        else
            iteration = DiagonalImplicit()
        end
    else
        iteration = FullImplicit()
    end
    return iteration
end

function fsal_prop(butcher::Matrix{T}) where T <: AbstractFloat
    ncol = size(butcher, 2)

    # TODO: if using implicit routine TRBDF2, then also need to check
    #       whether first stage is explicit in order to use f = f_tmp
    #       should probably test it out on TRBDF2 (fixed time step and embedded case)

    # remove any embedded rows
    B_square = butcher[1:ncol, :]

    # check if last, second-last rows are identical
    if B_square[end, :] == B_square[end-1, :]
        fsal = FSAL()
    else
        fsal = NotFSAL()
    end
    return fsal
end

function order_prop(name::Symbol, precision::Type{T}) where T <: AbstractFloat
    order = filter.(isdigit, split(string(name), "_"))
    order = parse.(precision, filter(x -> x != "", order))
    return SVector{length(order)}(order)
end

function explicit_stage_prop(butcher::Matrix{T}) where T <: AbstractFloat
    ncol = size(butcher, 2)
    A = butcher[1:(ncol-1), 2:end]
    stages = ncol-1
    explicit_stage = Vector{Bool}(undef, stages)
    for i in 1:stages
        explicit_stage[i] = all(x -> x == 0.0, A[i,i:end])
    end
    return explicit_stage |> SVector{stages}
end