"""
    debug_table(method::RungeKutta; tol_factor::Float64 = 1.7) 

Checks that the Butcher table satisfies the order conditions B[i,1] = ∑_{j>1} B[i,j]
within floating precision. By default, rows that fail the test are marked as `Broken`
(`Fail` if the error exceeds the tolerance by an order of magnitude). Primary and 
embedded iteration rows that fail the test are marked as `Broken` (`Fail` if the 
error exceeds the tolerance by the factor `tol_factor`).

Note: `tol_factor = 1.7` is the lowest factor I currently can reach
"""
function debug_table(method::RungeKutta; tol_factor::Float64 = 1.7)
    B = method.butcher

    for i in 1:size(B, 1)
        err = abs(B[i, 1] - sum(B[i, 2:end]))
        tol = norm(eps.(B[i, 2:end]), 1)        # sum floating precision errors
      
        if err > tol
            err = round(err |> Float64, sigdigits = 3)
            tol = round(tol |> Float64, sigdigits = 3)
            msg = "|B[$i,1] - ∑_{j>1} B[$i,j]| = $err > $tol. Check row $i in $method."

            # primary/embedded rows have a stricter test than stage rows
            if (err > tol_factor*tol && i >= size(B, 2)) || err > 10tol
                @error msg          # row should be fixed
                @test false         # test fail
            else
                @warn msg           # row might need fixing
                @test_broken false  # test broken
            end
        else
            @test true          # test passed
        end
    end
    nothing
end

function debug_iteration(method::RungeKutta, iteration::Type{<:Iteration})
    # TODO: test macro returns a message if failed
    prop = getproperty(method, :iteration)

    if prop isa iteration
        @test true 
    else 
        @error "$(method.name) iteration = $(typeof(prop)) but should be $iteration"
        @test false
    end
    nothing
end