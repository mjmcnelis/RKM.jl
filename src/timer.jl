"""
Sets a timer for the ODE solver.
"""
struct TimeLimit
    """Wall time (minutes)"""
    wtime_min::Float64
    """System time since epoch (seconds)"""
    time_sys::MVector{1,Float64}
    """Total number of time steps taken by the solver"""
    total_steps::MVector{1,Int64}
end

"""
    function TimeLimit(; wtime_min::Real = 60)

Outer constructor for `TimeLimit`.
"""
function TimeLimit(; wtime_min::Real = 60)
    @assert wtime_min >= 0 "wtime_min = $wtime_min is not greater than or equal to 0"
    time_sys = MVector{1,Float64}(time())
    total_steps = MVector{1,Int64}(0)

    return TimeLimit(wtime_min, time_sys, total_steps)
end

"""
    reset_timer!(timer::TimeLimit)

Resets the `timer` fields to the values given by the `TimeLimit` outer constructor

Required parameters: `timer`
"""
function reset_timer!(timer::TimeLimit)
    timer.time_sys[1] = time()
    timer.total_steps[1] = 0
    return nothing
end

"""
    continue_solver(t::Vector{T}, tf::T, timer::TimeLimit) where T <: AbstractFloat

Checks whether to continue running the ODE solver. The solver stops (`false`) if either
the simulation finishes `t >= tf` or the runtime exceeds the time limit set by `timer`.

Required parameters: `t`, `tf`, `timer`
"""
function continue_solver(t::Vector{T}, tf::T, timer::TimeLimit) where T <: AbstractFloat
    @unpack wtime_min, time_sys, total_steps = timer

    # note: check timer every 10 time steps
    if !isinf(wtime_min) && total_steps[1] % 10 == 0
        if time() > time_sys[1] + 60*wtime_min
            println("")
            # note: @code_warntype shows logger and err are type-unstable (from @warn)
            #       logger::Union{Nothing, Base.CoreLogging.AbstractLogger} (red)
            #       err::Any (red)
            # doesn't seem to hurt performance, would rather warn than print
            @warn "Exceeded time limit of $wtime_min minutes (stopping evolve_ode!...)\n"
            return false
        end
    end
    return t[1] < tf
end
