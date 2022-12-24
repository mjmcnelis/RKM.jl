using Revise
using RKM 
using Test
# TODO: export csv files for larger tables for viewing

floats = [Float32, Float64, BigFloat]

for precision in floats
    local RK_tables = vcat(
        get_runge_kutta_explicit_tables(; precision)...,
        get_runge_kutta_full_implicit_tables(; precision)...,
        get_runge_kutta_diagonal_implicit_tables(; precision)...
    )
    debug_table.(RK_tables)
    println("----------------------------------------------------------------------------")
end

println("\ndone\n")
