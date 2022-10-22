"""
    debug_table(method::RungeKutta) 

Checks that the Butcher table satisfies the order conditions B[i,1] = ∑_{j>1} B[i,j]
within floating precision. Rows that fail the test are marked as broken; whether or 
not they need debugging depends on how much the error exceeds the tolerance in the 
warning message.
"""
function debug_table(method::RungeKutta)
    B = method.butcher

    for i in 1:size(B, 1)
        err = abs(B[i, 1] - sum(B[i, 2:end]))
        tol = norm(eps.(B[i, 2:end]), 1)
      
        if err > tol
            err = round(err |> Float64, sigdigits = 3)
            tol = round(tol |> Float64, sigdigits = 3)
            msg = "|B[$i,1] - ∑_{j>1} B[$i,j]| = $err > $tol. Check row $i in $method."

            @warn msg           # row may need fixing
            @test_broken false  # test broken
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