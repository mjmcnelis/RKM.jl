
function method_is_fsal(butcher::SMatrix{N, M, T, NM}) where {N, M, T <: AbstractFloat, NM}
    ncol = size(butcher, 2)

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
