"""
Specifies the time evolution interval and initial time step.
"""
@kwdef struct TimeRange
    """Initial time"""
    t0::Float64
    """Final time"""
    tf::Float64
end

"""
Sets a timer for the ODE solver.
"""
struct TimeLimit
    """Wall time (minutes)"""
    wtime_min::Float64
    """System time since epoch (seconds)"""
    time_sys::MVector{1,Float64}
    """Counter for the number of time steps done so far by the solver"""
    counter::MVector{1,Int64}
end

"""
    function TimeLimit(; wtime_min::Real = 60)

Outer constructor for `TimeLimit`.
"""
function TimeLimit(; wtime_min::Real = 60)
    @assert wtime_min >= 0 "wtime_min = $wtime_min is not greater than or equal to 0"
    time_sys = MVector{1,Float64}(time())
    counter = MVector{1,Int64}(0)

    return TimeLimit(wtime_min, time_sys, counter)
end

"""
    reset_timer!(timer::TimeLimit)

Resets the `timer` fields to the values given by the `TimeLimit` outer constructor

Required parameters: `timer`
"""
function reset_timer!(timer::TimeLimit)
    timer.time_sys[1] = time()
    timer.counter[1] = 0
    return nothing
end

"""
    continue_solver(t::VectorMVector{1,T}, tf::T,
                    timer::TimeLimit) where T <: AbstractFloat

Checks whether to continue running the ODE solver. The solver stops (`false`) if either
the simulation finishes `t >= tf` or the runtime exceeds the time limit set by `timer`.

Required parameters: `t`, `tf`, `timer`
"""
function continue_solver(t::VectorMVector{1,T}, tf::T,
                         timer::TimeLimit) where T <: AbstractFloat
    @unpack counter, wtime_min, time_sys = timer

    counter[1] += 1
    # note: check timer every 10 time steps
    if !isinf(wtime_min) && counter[1] % 10 == 0
        if time() > time_sys[1] + 60*wtime_min
            println("")
            @warn "Exceeded time limit of $wtime_min minutes (stopping evolve_ode!...)\n"
            return false
        end
    end
    return t[1] < tf
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
