
abstract type Wrapper end

struct ODEWrapper{T, P, F} <: Wrapper where {T <: AbstractFloat, F <: Function}
    t::Vector{T}
    p::Vector{Float64}
    abstract_params::P
    dy_dt!::F
end

# used by ForwardDiff: ode_wrap!(f, y)
function (ode_wrap!::ODEWrapper)(f, y)
    @unpack t, p, abstract_params = ode_wrap!
    ode_wrap!.dy_dt!(f, y, t[1]; p, abstract_params)
end

# used by RKM routines: ode_wrap!(f, t, y)
function (ode_wrap!::ODEWrapper)(f, t, y)
    @unpack p, abstract_params = ode_wrap!
    ode_wrap!.dy_dt!(f, y, t; p, abstract_params)
end
