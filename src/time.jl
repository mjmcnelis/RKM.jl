"""
Specifies the time evolution interval and initial time step.
"""
@kwdef struct TimeRange
    """Initial time"""
    t0::Float64
    """Final time"""
    tf::Float64
    """Initial time step"""
    dt0::Float64
end

"""
Sets a timer for the ODE solver.
"""
struct TimeLimit
    """Wall time time (minutes)"""
    wtime_min::Int64
    """Time limit in terms of DateTime"""
    time_limit::DateTime
    """Check timer after this number of time steps"""
    frequency::Int64
    """Counter for number of time steps done so far by the solver"""
    counter::MVector{1,Int64}
end

"""
    function TimeLimit(; wtime_min::Int64 = 60, frequency::Int64 = 100)

Outer constructor for `TimeLimit`.
"""
function TimeLimit(; wtime_min::Int64 = 60, frequency::Int64 = 100)
    time_limit = now() + Minute(round(wtime_min))
    counter = MVector{1,Int64}(0)

    TimeLimit(wtime_min, time_limit, frequency, counter)
end

"""
    continue_solver(t::MVector{1,T}, tf::T2, 
                    timer::TimeLimit) where {T <: AbstractFloat, T2 <: AbstractFloat}

Checks whether to continue running the ODE solver. The solver stops (`false`) if either
the simulation finishes `t >= tf` or the runtime exceeds the time limit set by `timer`.

Required parameters: `t`, `tf`, `timer`
"""
function continue_solver(t::MVector{1,T}, tf::T2, 
                         timer::TimeLimit) where {T <: AbstractFloat, T2 <: AbstractFloat}
    t[1] < tf && !past_time_limit(timer)
end

"""
    past_time_limit(timer::TimeLimit)

Checks whether the runtime has exceeded the time limit set by `timer`.

Required parameters: `timer`
"""
function past_time_limit(timer::TimeLimit)
    @unpack counter, frequency, time_limit = timer

    @.. counter += 1
    if counter[1] % frequency == 0 && now() > time_limit
        @warn "\nExceeded time limit of $(timer.wtime_min) minutes (stop evolve...)\n"
        true
    else
        false
    end
end
