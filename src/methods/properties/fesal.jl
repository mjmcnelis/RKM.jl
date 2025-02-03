"""
    get_fesal(A_T::SMatrix{S, S, T, S2},
              b::SVector{S, T}) where {S, S2, T <: AbstractFloat}

Checks whether the Runge-Kutta method has the First Explicit Same As Last (FESAL) property.
The boolean `fesal` is set true if the last stage coefficients `A_T[:,end]` are identical
to the primary update coefficients `b`.

FESAL drops the requirement that the first stage is explicit, which makes it more relaxed
the conventional First Same As Last (FSAL) property. Some implicit Runge-Kutta methods
(e.g. BackwardEuler1) are allowed to reuse the last stage for Hermite interpolation or an
initial guess for the next implicit stage.

Required parameters: `A_T`, `b`
"""
function get_fesal(A_T::SMatrix{S, S, T, S2},
                   b::SVector{S, T}) where {S, S2, T <: AbstractFloat}
    fesal = all(A_T[:,end] .== b) ? true : false
    return fesal
end
