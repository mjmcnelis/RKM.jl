using Revise
using Test

@info "Starting precompile test..."

stats = @timed using RKM

println("\tPrecompile time = $(round(stats.time, digits = 2)) seconds")
@test stats.time < 10.0

@info "...done"
println("")
