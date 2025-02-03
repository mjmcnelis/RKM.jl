
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
