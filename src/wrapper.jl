
abstract type Wrapper end

struct ODEWrapperState{T, P, F} <: Wrapper where {T <: AbstractFloat, F <: Function}
    t::Vector{T}
    p::Vector{Float64}
    abstract_params::P
    dy_dt!::F
    FE::MVector{1,Int64}
    # JE::Vector{Int64}
end

function ODEWrapperState(t, p, abstract_params, dy_dt!)#, method)
    FE = MVector{1,Int64}(0)
    # TODO: not sure why I couldn't do MVector
    # JE = zeros(Int64, method.stages)
    return ODEWrapperState(t, p, abstract_params, dy_dt!, FE)#, JE)
end

# note: t is a vector in first method but float in the other
#       is that a type instability conflict?

# used by RKM routines: ode_wrap!(f, t, y)
function (ode_wrap!::ODEWrapperState)(f::Vector{T}, t::T,
                                      y::Vector{T}) where T <: AbstractFloat
    p = ode_wrap!.p
    abstract_params = ode_wrap!.abstract_params
    FE = ode_wrap!.FE

    FE[1] += 1
    ode_wrap!.dy_dt!(f, y, t; p, abstract_params)

    return nothing
end

# used by FiniteDiff and ForwardDiff: ode_wrap!(f, y)
function (ode_wrap!::ODEWrapperState)(f::Vector{R}, y::Vector{R}) where R <: Real
    t = ode_wrap!.t
    p = ode_wrap!.p
    abstract_params = ode_wrap!.abstract_params
    FE = ode_wrap!.FE

    FE[1] += 1
    ode_wrap!.dy_dt!(f, y, t[1]; p, abstract_params)

    return nothing
end

struct ODEWrapperParam{T, P, F} <: Wrapper where {T <: AbstractFloat, F <: Function}
    t::Vector{T}
    y::Vector{T}    # new vector or reuse one?
    abstract_params::P
    dy_dt!::F
    FE::MVector{1,Int64}
    # JE::Vector{Int64}
end

function ODEWrapperParam(t, y, abstract_params, dy_dt!)#, method)
    FE = MVector{1,Int64}(0)
    # JE = zeros(Int64, method.stages)
    return ODEWrapperParam(t, y, abstract_params, dy_dt!, FE)#, JE)
end

# used by FiniteDiff and ForwardDiff: ode_wrap!(f, p)
# TODO: does not work for R = DoubleFloat unless convert p
function (ode_wrap!::ODEWrapperParam)(f::Vector{R}, p::Vector{R}) where R <: Real
    y = ode_wrap!.y
    t = ode_wrap!.t
    abstract_params = ode_wrap!.abstract_params
    FE = ode_wrap!.FE

    FE[1] += 1
    ode_wrap!.dy_dt!(f, y, t[1]; p, abstract_params)

    return nothing
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