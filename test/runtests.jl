using SafeTestsets
using Test

@testset "RKM.jl tests" begin
    @safetestset "Order condition test" begin include("butcher_test.jl") end
    @safetestset "Iteration test" begin include("iteration_test.jl") end
end