"""
    debug_table(method::RungeKutta; p::Int = 1) 

Checks that the butcher table satisfies the order conditions B[i,1] = ∑_{j>1} B[i,j]
within floating precision. If the last row fails the test (strict), a failed test is 
registered. Other rows failing the test are marked as broken but may need debugging.
"""
function debug_table(method::RungeKutta; p::Int = 1)
    B = method.butcher

    for i in 1:size(B, 1)
        err = abs(B[i, 1] - sum(B[i, 2:end]))
        tol = norm(eps.(B[i, 2:end]), p)        # p = 1 seems to work best
      
        if err > tol
            err = round(err |> Float64, sigdigits = 3)
            tol = round(tol |> Float64, sigdigits = 3)
            msg = "|B[$i,1] - ∑_{j>1} B[$i,j]| = $err > $tol. Fix row $i in $method."

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