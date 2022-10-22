
# TODO: make docstring
function get_all_runge_kutta_tables(; precision::Type{<:AbstractFloat})
    vcat(
        get_runge_kutta_explicit_tables(; precision)...,
        get_runge_kutta_full_implicit_tables(; precision)...,
        get_runge_kutta_diagonal_implicit_tables(; precision)...,
    )
end