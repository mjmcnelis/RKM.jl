using Revise
using RKM
using DoubleFloats
using Test
include(joinpath(RKM_root, "test/tables/get_explicit_tables.jl"))
include(joinpath(RKM_root, "test/tables/get_implicit_tables.jl"))
include(joinpath(RKM_root, "test/utils.jl"))

floats = [Float32, Float64, Double64, BigFloat]

for precision in floats
    local RK_explicit      = get_runge_kutta_explicit_tables(; precision)
    local RK_full_implicit = get_runge_kutta_full_implicit_tables(; precision)
    local RK_diag_implicit = get_runge_kutta_diagonal_implicit_tables(; precision)

    debug_iteration.(RK_explicit,      Explicit)
    debug_iteration.(RK_full_implicit, FullImplicit)
    debug_iteration.(RK_diag_implicit, DiagonalImplicit)
end

println("\ndone\n")
