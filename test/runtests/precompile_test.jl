using Test

stats = @timed using RKM 
@show stats.time 

# 3/28/23: current time is 2.7s (2.0s if remove tables from src/)

# want precompile time to be shorter than 10s
@test stats.time < 10.0

println("\ndone")
