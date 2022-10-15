
# TODO: break up properties in directory 

abstract type Iteration end 
abstract type Implicit <: Iteration end

struct Explicit <: Iteration end
struct DiagonalImplicit <: Implicit end
struct FullImplicit <: Implicit end

abstract type FirstSameAsLast end 
struct FSAL <: FirstSameAsLast end 
struct NoFSAL <: FirstSameAsLast end

# TODO: determine if butcher table fixed (square matrix) or embedded (not square matrix)

function iteration_prop(butcher::Matrix{<:AbstractFloat})
    ncol = size(butcher, 2)
    A = butcher[1:(ncol-1), 2:end]      # get sub-matrix A_{ij}

    if tril(A) == A                     # check if A is lower triangular
        if all(x -> x == 0, diag(A))    # check if diagonal elements are all zero
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

    fsal = FSAL()

    fsal
end