
abstract type Wrapper end

struct ODEWrapperState{T, P, F} <: Wrapper where {T <: AbstractFloat, F <: Function}
    t::Vector{T}
    p::Vector{Float64}
    abstract_params::P
    dy_dt!::F
end

# note: t is a vector in first method but float in the other
#       is that a type instability conflict?

# used by RKM routines: ode_wrap!(f, t, y)
function (ode_wrap!::ODEWrapperState)(f::Vector{T}, t::T,
                                      y::Vector{T}) where T <: AbstractFloat
    @unpack p, abstract_params = ode_wrap!
    ode_wrap!.dy_dt!(f, y, t; p, abstract_params)
end

# used by FiniteDiff and ForwardDiff: ode_wrap!(f, y)
function (ode_wrap!::ODEWrapperState)(f::Vector{R}, y::Vector{R}) where R <: Real
    @unpack t, p, abstract_params = ode_wrap!
    ode_wrap!.dy_dt!(f, y, t[1]; p, abstract_params)
end

struct ODEWrapperParam{T, P, F} <: Wrapper where {T <: AbstractFloat, F <: Function}
    t::Vector{T}
    y::Vector{T}    # new vector or reuse one?
    abstract_params::P
    dy_dt!::F
end

# used by FiniteDiff and ForwardDiff: ode_wrap!(f, p)
function (ode_wrap!::ODEWrapperParam)(f::Vector{R}, p::Vector{R}) where R <: Real
    @unpack y, t, abstract_params = ode_wrap!
    ode_wrap!.dy_dt!(f, y, t[1]; p, abstract_params)
end
