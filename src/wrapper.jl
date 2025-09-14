
abstract type Wrapper end

struct ODEWrapperState{T, P, F} <: Wrapper where {T <: AbstractFloat, F <: Function}
    t::Vector{T}
    p::Vector{T}
    abstract_params::P
    dy_dt!::F
    evaluations::MVector{2,Int64}
    time_subroutine::Bool
    runtime::MVector{1,Float64}
end

function ODEWrapperState(t, p, abstract_params, dy_dt!, time_subroutine)
    evaluations = MVector{2,Int64}(0, 0)
    runtime = MVector{1,Float64}(0.0)

    return ODEWrapperState(t, p, abstract_params, dy_dt!, evaluations,
                           time_subroutine, runtime)
end

# note: t is a vector in first method but float in the other
#       is that a type instability conflict?

# used by RKM routines: ode_wrap!(f, t, y)
function (ode_wrap!::ODEWrapperState)(f::Vector{T}, t::T,
                                      y::Vector{T}) where T <: AbstractFloat
    p = ode_wrap!.p
    abstract_params = ode_wrap!.abstract_params
    evaluations = ode_wrap!.evaluations
    time_subroutine = ode_wrap!.time_subroutine
    runtime = ode_wrap!.runtime

    if time_subroutine && evaluations[1] % SAMPLE_INTERVAL == 0
        stats = @timed ode_wrap!.dy_dt!(f, y, t; p, abstract_params)
        runtime[1] += SAMPLE_INTERVAL*stats.time
    else
        ode_wrap!.dy_dt!(f, y, t; p, abstract_params)
    end
    evaluations[1] += 1

    return nothing
end

# used by FiniteDiff and ForwardDiff: ode_wrap!(f, y)
function (ode_wrap!::ODEWrapperState)(f::Vector{R}, y::Vector{R}) where R <: Real
    t = ode_wrap!.t
    p = ode_wrap!.p
    abstract_params = ode_wrap!.abstract_params
    evaluations = ode_wrap!.evaluations

    ode_wrap!.dy_dt!(f, y, t[1]; p, abstract_params)
    evaluations[2] += 1

    return nothing
end

struct ODEWrapperParam{T, P, F} <: Wrapper where {T <: AbstractFloat, F <: Function}
    t::Vector{T}
    y::Vector{T}    # new vector or reuse one?
    abstract_params::P
    dy_dt!::F
    evaluations::MVector{1,Int64}
end

function ODEWrapperParam(t, y, abstract_params, dy_dt!)
    evaluations = MVector{1,Int64}(0)

    return ODEWrapperParam(t, y, abstract_params, dy_dt!, evaluations)
end

# used by FiniteDiff and ForwardDiff: ode_wrap!(f, p)
function (ode_wrap!::ODEWrapperParam)(f::Vector{R}, p::Vector{R}) where R <: Real
    y = ode_wrap!.y
    t = ode_wrap!.t
    abstract_params = ode_wrap!.abstract_params
    evaluations = ode_wrap!.evaluations

    ode_wrap!.dy_dt!(f, y, t[1]; p, abstract_params)
    evaluations[1] += 1

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