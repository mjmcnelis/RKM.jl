
abstract type Wrapper end

struct ODEWrapper{T, P, F} <: Wrapper where {T <: AbstractFloat, F <: Function}
    t::VectorMVector{1,T}
    p::P
    dy_dt!::F
end

(ode_wrap!::ODEWrapper)(f, y) = ode_wrap!.dy_dt!(f, y; t = ode_wrap!.t[1], p = ode_wrap!.p)
(ode_wrap!::ODEWrapper)(f, t, y) = ode_wrap!.dy_dt!(f, y; t, p = ode_wrap!.p)
