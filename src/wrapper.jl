
abstract type Wrapper end

# TODO: add parameter list
struct ODEWrapper{T, F} <: Wrapper where {T <: AbstractFloat, F <: Function}
    t::VectorMVector{1,T}
    dy_dt!::F
end

(ode_wrap::ODEWrapper)(f, y) = ode_wrap.dy_dt!(f, ode_wrap.t[1], y)