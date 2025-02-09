"""
    get_fesal(A_T::SMatrix{S, S, T, S2}, b::SVector{S, T},
              c::SVector{S, T}) where {S, S2, T <: AbstractFloat}

Checks whether the Runge-Kutta method has the First Explicit Same As Last (FESAL) property.
The boolean `fesal` is set true if the last stage is equal to an explicit first stage
evaluated at the next time step.

Unlike the conventional First Same As Last (FSAL) property, FESAL does not require that
the first stage of the Butcher tableau is explicit. Therefore, some implicit Runge-Kutta
methods like BackwardEuler1 are allowed to reuse the last stage for Hermite interpolation.

Required parameters: `A_T`, `b`, `c`
"""
function get_fesal(A_T::SMatrix{S, S, T, S2}, b::SVector{S, T},
                   c::SVector{S, T}) where {S, S2, T <: AbstractFloat}
    fesal = all(A_T[:,end] .== b) && c[end] == 1.0 ? true : false
    return fesal
end
