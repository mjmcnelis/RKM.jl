"""
    debug_table(method::RungeKutta) 

Checks that the butcher table satisfies the order conditions B[i,1] = ∑_{j>1} B[i,j]
within floating precision. If the last row fails the test (strict), a failed test is 
registered. Other rows failing the test are marked as broken but may need debugging.
"""
function debug_table(method::RungeKutta)
    B = method.butcher
    for i in 1:size(B, 1)
        err = B[i, 1] - sum(B[i, 2:end])
        # TODO: adjust tolerance depending on Float type
        if abs(err) > 1e-15
            msg = "B[$i,1] - ∑_{j>1} B[$i,j] = $err. Fix row $i in $method."
            if i == size(B, 1)
                @error msg          # row must be fixed
                @test false         # test failed
            else
                @warn msg           # row may need fixing
                @test_broken false  # test broken
            end 
        else
            @test true              # test passed
        end
    end
    nothing
end
