
abstract type Wrapper end

struct ODEWrapperState{T, P, F} <: Wrapper where {T <: AbstractFloat, F <: Function}
    t::Vector{T}
    p::Vector{Float64}
    abstract_params::P
    dy_dt!::F
    FE::MVector{1,Int64}
end

function ODEWrapperState(t, p, abstract_params, dy_dt!)
    FE = MVector{1,Int64}(0)
    return ODEWrapperState(t, p, abstract_params, dy_dt!, FE)
end

# note: t is a vector in first method but float in the other
#       is that a type instability conflict?

# used by RKM routines: ode_wrap!(f, t, y)
function (ode_wrap!::ODEWrapperState)(f::Vector{T}, t::T,
                                      y::Vector{T}) where T <: AbstractFloat
    @unpack p, abstract_params, FE = ode_wrap!
    FE[1] += 1
    ode_wrap!.dy_dt!(f, y, t; p, abstract_params)
end

# used by FiniteDiff and ForwardDiff: ode_wrap!(f, y)
function (ode_wrap!::ODEWrapperState)(f::Vector{R}, y::Vector{R}) where R <: Real
    @unpack t, p, abstract_params, FE = ode_wrap!
    FE[1] += 1
    ode_wrap!.dy_dt!(f, y, t[1]; p, abstract_params)
end

struct ODEWrapperParam{T, P, F} <: Wrapper where {T <: AbstractFloat, F <: Function}
    t::Vector{T}
    y::Vector{T}    # new vector or reuse one?
    abstract_params::P
    dy_dt!::F
    FE::MVector{1,Int64}
end

function ODEWrapperParam(t, y, abstract_params, dy_dt!)
    FE = MVector{1,Int64}(0)
    return ODEWrapperParam(t, y, abstract_params, dy_dt!, FE)
end

# used by FiniteDiff and ForwardDiff: ode_wrap!(f, p)
# TODO: does not work for R = DoubleFloat unless convert p
function (ode_wrap!::ODEWrapperParam)(f::Vector{R}, p::Vector{R}) where R <: Real
    @unpack y, t, abstract_params, FE = ode_wrap!
    FE[1] += 1
    ode_wrap!.dy_dt!(f, y, t[1]; p, abstract_params)
end

function set_wrapper!(ode_wrap!::ODEWrapperState, t::T) where T <: AbstractFloat
    ode_wrap!.t[1] = t
    return nothing
end

function set_wrapper!(ode_wrap!::ODEWrapperParam,
                      t::T, y::Vector{T}) where T <: AbstractFloat
    ode_wrap!.t[1] = t
    # note: safe but may be redundant if they share same pointer,
    # depends on what I need to set ode_wrap!.y to (so far y_tmp)
    @.. ode_wrap!.y = y
    return nothing
end