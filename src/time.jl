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
    reset_timer(timer::TimeLimit)

Constructs a new `TimeLimit` from `timer` by updating `time_limit`
and `counter` fields.

Required parameters: `timer`
"""
function reset_timer(timer::TimeLimit)
    @unpack wtime_min, frequency = timer 
    return TimeLimit(; wtime_min, frequency)
end

"""
    monitor_progess(t::Union{Vector{T}, MVector{1,T}}, progress::Progress,
                    checkpoints::Vector{T}) where T <: AbstractFloat

Updates the `progress` meter percentage points depending on how many 
`checkpoints` the current time `t` has passed.

Required parameters: `t, `progress`, `checkpoints`
"""
function monitor_progess(t::Union{Vector{T}, MVector{1,T}}, progress::Progress,
                         checkpoints::Vector{T}) where T <: AbstractFloat
    if length(checkpoints) > 1
        dt = checkpoints[2] - checkpoints[1]
        idx = Int(floor(Float64((t[1] - checkpoints[1])/dt))) + 1
        idx = min(length(checkpoints), idx)
        for i = 1:idx
            next!(progress)
            popfirst!(checkpoints)
        end  
    else 
        if t[1] >= checkpoints[1]
            next!(progress)
            sleep(1e-3)
        end
    end
    return nothing
end

"""
    continue_solver(t::Union{Vector{T}, MVector{1,T}}, tf::T,
                    timer::TimeLimit) where T <: AbstractFloat

Checks whether to continue running the ODE solver. The solver stops (`false`) if either
the simulation finishes `t >= tf` or the runtime exceeds the time limit set by `timer`.

Required parameters: `t`, `tf`, `timer`
"""
function continue_solver(t::Union{Vector{T}, MVector{1,T}}, tf::T,
                         timer::TimeLimit) where T <: AbstractFloat
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
    # check timer for every N = frequency time steps
    if counter[1] % frequency == 0
        if now() > time_limit
            @warn "\nExceeded time limit of $(timer.wtime_min) minutes (stop evolve...)\n"
            return true
        end
    end 
    return false
end
