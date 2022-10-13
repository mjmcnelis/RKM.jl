using SafeTestsets
using Test

@testset "RKM.jl tests" begin
    @safetestset "Butcher tables" begin include("butcher_test.jl") end
end