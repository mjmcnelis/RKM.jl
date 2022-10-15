
# TODO: break up properties in directory 

abstract type Iteration end 
abstract type Implicit <: Iteration end

struct Explicit <: Iteration end
struct DiagonalImplicit <: Implicit end
struct FullImplicit <: Implicit end

abstract type FirstSameAsLast end 
struct FSAL <: FirstSameAsLast end 
struct NotFSAL <: FirstSameAsLast end

# TODO: determine if butcher table fixed (square matrix) or embedded (not square matrix)

function iteration_prop(butcher::Matrix{<:AbstractFloat})
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
    iteration
end

function fsal_prop(butcher::Matrix{<:AbstractFloat})
    ncol = size(butcher, 2)

    # remove any embedded rows
    B_square = butcher[1:ncol, :]

    # check if last, second-last rows are identical
    if B_square[end, :] == B_square[end-1, :]
        fsal = FSAL() 
    else
        fsal = NotFSAL()
    end
    fsal
end