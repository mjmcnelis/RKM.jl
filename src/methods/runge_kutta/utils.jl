
function debug_table(method::RungeKutta)
    B = method.butcher
    # check B[i,1] = ∑_{j>1} B[i,j] (strict condition for last row)
    for i in 1:size(B, 1)
        err = B[i, 1] - sum(B[i, 2:end])
        # TODO: adjust tolerance depending Float type
        if abs(err) > 1.0e-15
            @warn "B_i - ∑_{j>1} B_ij = $err. Fix row $i in $method."
            if i == size(B, 1)
                # test fails (strict condition for last row)
                return false
            end
        end
    end
    # test passed
    true
end