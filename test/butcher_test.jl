using Revise
using RKM
using DoubleFloats
using Test
include(joinpath(RKM_root, "test/tables/get_explicit_tables.jl"))
include(joinpath(RKM_root, "test/tables/get_implicit_tables.jl"))
include(joinpath(RKM_root, "test/utils.jl"))
# TODO: export csv files for larger tables for viewing

for precision in [Float32, Float64, Double64, BigFloat]
    RK_tables = vcat(
        get_runge_kutta_explicit_tables(; precision)...,
        get_runge_kutta_full_implicit_tables(; precision)...,
        get_runge_kutta_diagonal_implicit_tables(; precision)...
    )
    debug_table.(RK_tables)
    println("----------------------------------------------------------------------------")
end

println("\ndone\n")
