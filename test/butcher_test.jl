using Revise
using RKM 
using Test
# TODO: export csv files for larger tables for viewing

floats = [Float32, Float64, BigFloat]

for precision in floats
    local RK_tables = get_all_runge_kutta_tables(; precision)
    debug_table.(RK_tables)
    println("----------------------------------------------------------------------------")
end

println("\ndone\n")
