using Test

stats = @timed using RKM 
@show stats.time 

# want precompile time to be shorter than 10s
@test stats.time < 10.0

println("\ndone")
