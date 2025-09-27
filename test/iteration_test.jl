using RKM
using DoubleFloats
using Test
include(joinpath(RKM_root, "test/runtests/tables/get_explicit_tables.jl"))
include(joinpath(RKM_root, "test/runtests/tables/get_implicit_tables.jl"))
include(joinpath(RKM_root, "test/runtests/utils.jl"))

RK_explicit      = get_runge_kutta_explicit_tables(; precision = Float64)
RK_full_implicit = get_runge_kutta_full_implicit_tables(; precision = Float64)
RK_diag_implicit = get_runge_kutta_diagonal_implicit_tables(; precision = Float64)

println("\ndone\n")
