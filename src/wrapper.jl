
abstract type Wrapper end

struct ODEWrapper{T, P, F} <: Wrapper where {T <: AbstractFloat, F <: Function}
    t::VectorMVector{2,T}   # note: for now, 2nd index is not being used/updated
    p::P
    dy_dt!::F
end

# used by ForwardDiff: ode_wrap!(f, y)
(ode_wrap!::ODEWrapper)(f, y) = ode_wrap!.dy_dt!(f, y; t = ode_wrap!.t[1], p = ode_wrap!.p)

# used by RKM routines: ode_wrap!(f, t, y)
(ode_wrap!::ODEWrapper)(f, t, y) = ode_wrap!.dy_dt!(f, y; t, p = ode_wrap!.p)
