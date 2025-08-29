
abstract type Wrapper end

struct ODEWrapperState{T, P, F} <: Wrapper where {T <: AbstractFloat, F <: Function}
    t::Vector{T}
    p::Vector{Float64}
    abstract_params::P
    dy_dt!::F
    evaluations::MVector{1,Int64}
    benchmarks::Bool
    subroutine_time::MVector{1,Float64}
end

function ODEWrapperState(t, p, abstract_params, dy_dt!, benchmarks)
    evaluations = MVector{1,Int64}(0)
    subroutine_time = MVector{1,Float64}(0.0)

    return ODEWrapperState(t, p, abstract_params, dy_dt!, evaluations,
                           benchmarks, subroutine_time)
end

# note: t is a vector in first method but float in the other
#       is that a type instability conflict?

# used by RKM routines: ode_wrap!(f, t, y)
function (ode_wrap!::ODEWrapperState)(f::Vector{T}, t::T,
                                      y::Vector{T}) where T <: AbstractFloat
    p = ode_wrap!.p
    abstract_params = ode_wrap!.abstract_params
    evaluations = ode_wrap!.evaluations
    benchmarks = ode_wrap!.benchmarks
    subroutine_time = ode_wrap!.subroutine_time

    if benchmarks && evaluations[1] % SAMPLE_INTERVAL == 0
        stats = @timed ode_wrap!.dy_dt!(f, y, t; p, abstract_params)
        subroutine_time[1] += SAMPLE_INTERVAL*stats.time
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

    # could do evaluations_J[1] += 1 here
    # TMP comment out b/c it interferes with subroutine_time accumulation
    # evaluations[1] += 1
    ode_wrap!.dy_dt!(f, y, t[1]; p, abstract_params)

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
# TODO: does not work for R = DoubleFloat unless convert p
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