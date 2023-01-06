import LinearAlgebra: norm

"""
    debug_table(method::RungeKutta; tol_fact_iter = 1.86, tol_fact_stage = 10.0)

Checks that the Butcher table satisfies the order conditions B[i,1] = ∑_{j>1} B[i,j]
within floating precision. Stage rows that fail the test are marked as `Broken`
(`Fail` if the error exceeds the tolerance by a factor of `tol_fact_stage`). Primary
and embedded iteration rows that fail the test are marked as `Broken` (`Fail` if the
error exceeds the tolerance by the factor `tol_fact_iter`).

Note: `tol_fact_iter = 1.86` is the lowest factor I currently can reach.
"""
function debug_table(method::RungeKutta; tol_fact_iter = 1.86, tol_fact_stage = 10.0)
    B = method.butcher

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

function debug_iteration(method::RungeKutta, iteration::Type{I}) where I <: Iteration
    # TODO: test macro returns a message if failed
    prop = getproperty(method, :iteration)

    if !(prop isa iteration)
        @error "$(method.name) iteration = $(typeof(prop)) but should be $iteration"
        @test false
    end
    nothing
end
