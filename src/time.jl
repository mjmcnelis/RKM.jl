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
    """Wall time (minutes)"""
    wtime_min::Int64
    """Determines whether or not the runtime has exceeded the time limit"""
    past_time::MVector{1,Bool}
    """Determines whether or not the ODE solver has finished running""" 
    solver_finished::MVector{1,Bool}
    """Counter for the number of time steps done so far by the solver"""
    counter::MVector{1,Int64}
end

"""
    function TimeLimit(; wtime_min::Int64 = 60)

Outer constructor for `TimeLimit`.
"""
function TimeLimit(; wtime_min::Int64 = 60)
    past_time = MVector{1,Bool}(false)
    solver_finished = MVector{1,Bool}(false)
    counter = MVector{1,Int64}(0)

    return TimeLimit(wtime_min, past_time, solver_finished, counter)
end

"""
    reset_timer!(timer::TimeLimit)

Resets the `timer` fields to values given by `TimeLimit` outer constructor

Required parameters: `timer`
"""
function reset_timer!(timer::TimeLimit)
    timer.past_time[1] = false 
    timer.solver_finished[1] = false
    timer.counter[1] = 0
    return nothing
end

"""
    start_timer!(timer::TimeLimit)

Starts the `timer` and set the field `past_time` to 
`true` when the runtime exceeds the time limit.

Required parameters: `timer`
"""
function start_timer!(timer::TimeLimit)
    @unpack wtime_min, past_time, solver_finished = timer 
    @async begin 
        # @info "Timer set to $wtime_min minute(s)"
        sleep(60 * wtime_min)
        if !solver_finished[1]
            @warn "\nExceeded time limit of $wtime_min minutes (stopping evolve_ode...)\n"
        end
        past_time[1] = true
    end 
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
    # TODO: add counter somewhere else? 
    timer.counter[1] += 1
    if t[1] >= tf 
        timer.solver_finished[1] = true 
    end
    return t[1] < tf && !timer.past_time[1]
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
