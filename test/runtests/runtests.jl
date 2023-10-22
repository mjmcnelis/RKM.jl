using SafeTestsets
using Test

@testset "RKM.jl tests" begin
    @safetestset "Precompile test" begin include("precompile_test.jl") end
    @safetestset "Order condition test" begin include("butcher_test.jl") end
    @safetestset "Iteration test" begin include("iteration_test.jl") end
    @safetestset "Allocation test" begin include("allocation_test.jl") end
end
