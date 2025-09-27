
function reconstruct_butcher(method::RungeKutta)

    stages = method.stages
    c = method.c
    A_T = method.A_T
    b = method.b
    b_hat = method.b_hat

    precision = typeof(c[1])

    ncol = stages + 1
    nrow = b_hat == b ? ncol : ncol + 1
    butcher = zeros(precision, nrow, ncol)

    butcher[1:ncol-1, 2:ncol] .= A_T'
    butcher[1:ncol-1, 1] .= c
    butcher[ncol, 2:ncol] .= b
    butcher[nrow, 2:ncol] .= b_hat
    butcher[ncol, 1] = 1.0
    butcher[nrow, 1] = 1.0

    return butcher
end

"""
    debug_table(method::RungeKutta; tol_fact_iter = 1.86, tol_fact_stage = 10.0)

Checks that the Butcher table satisfies the order conditions B[i,1] = ∑ⱼ₌₁  B[i,j]
within floating precision. Stage rows that fail the test are marked as `Broken`
(`Fail` if the error exceeds the tolerance by a factor of `tol_fact_stage`). Primary
and embedded iteration rows that fail the test are marked as `Broken` (`Fail` if the
error exceeds the tolerance by the factor `tol_fact_iter`).

Note: `tol_fact_iter = 1.86` is the lowest factor we currently can reach.
"""
function debug_table(method::RungeKutta; tol_fact_iter = 1.86, tol_fact_stage = 10.0)
    B = reconstruct_butcher(method)

    for i in 1:size(B, 1)
        err = abs(B[i, 1] - sum(B[i, 2:end]))
        tol = norm(eps.(B[i, 2:end]), 1)        # sum floating precision errors

        if err > tol
            err = round(err |> Float64, sigdigits = 3)
            tol = round(tol |> Float64, sigdigits = 3)

            name = replace(String(method.name), "_" => "")
            msg = "|B[$i,1] - ∑_{j>1} B[$i,j]| = $err > $tol. Check row $i in $name."

            # primary/embedded iteration rows should have a stricter test than stage rows
            if (err > tol_fact_iter*tol && i >= size(B, 2)) || err > tol_fact_stage*tol
                @error msg          # row should be fixed
                @test false         # test fail
                # @test_broken false    # for debugging
            else
                @warn msg           # row might need fixing
                @test_broken false  # test broken
            end
        end
    end
    nothing
end
