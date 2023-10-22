using Revise
using Test

stats = @timed using RKM
@info "Precompile time = $(stats.time) seconds"

# 3/28/23: current time is 2.7s (2.0s if remove tables from src/)

# idealy want precompile time to be shorter than 10s
@test stats.time < 10.0

# after adding SciMLBase and LinearSolve for implicit routines, precompile time ~ 16s
# @test stats.time < 20.0

println("\ndone")
