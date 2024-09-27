using Test
@info "Starting precompile test..."

stats = @timed using RKM
println("\tPrecompile time = $(round(stats.time, digits = 2)) seconds")
if stats.time > 10.0
    @test_broken false
end
@test stats.time < 20.0

@info "...done"
println("")
