"""
Sets a timer for the ODE solver.
"""
struct TimeLimit
    """Wall time (minutes)"""
    wtime_min::Float64
    """Initial and current system time since epoch (seconds)"""
    time_sys::MVector{2,Float64}
    """Runtime in HH:MM:SS format"""
    runtime::Vector{String}
    """Runtime that was previously displayed on progress bar (seconds)"""
    runtime_prev::MVector{1,Int64}
    """Dynamic switch to update progress bar every second in real time"""
    display_values::MVector{1,Bool}
    """Total number of time steps taken by the solver"""
    total_steps::MVector{1,Int64}
end

"""
    function TimeLimit(; wtime_min::Real = 60)

Outer constructor for `TimeLimit`.
"""
function TimeLimit(; wtime_min::Real = 60)
    @assert wtime_min >= 0 "wtime_min = $wtime_min is not greater than or equal to 0"
    t_epoch = time()
    time_sys = MVector{2,Float64}(t_epoch, t_epoch)
    runtime = String["00:00:00"]
    runtime_prev = MVector{1,Int64}(0)
    display_values = MVector{1,Bool}(false)
    total_steps = MVector{1,Int64}(0)

    return TimeLimit(wtime_min, time_sys, runtime, runtime_prev,
                     display_values, total_steps)
end

"""
    reset_timer!(timer::TimeLimit)

Resets the `timer` fields to the values given by the `TimeLimit` outer constructor

Required parameters: `timer`
"""
function reset_timer!(timer::TimeLimit)
    t_epoch = time()
    timer.time_sys .= t_epoch
    timer.runtime[1] = "00:00:00"
    timer.runtime_prev[1] = 0
    timer.display_values[1] = false
    timer.total_steps[1] = 0
    return nothing
end

# TODO: make docstring
function set_current_system_time!(timer::TimeLimit)
    timer.time_sys[2] = time()
    return nothing
end

# TODO: make docstring
function set_runtime_display!(timer::TimeLimit)
    time_sys = timer.time_sys
    runtime = timer.runtime
    runtime_prev = timer.runtime_prev
    display_values = timer.display_values

    dt_sys = floor(Int64, time_sys[2] - time_sys[1])

    if dt_sys > runtime_prev[1]
        runtime_prev[1] = dt_sys
        h = floor(Int64, dt_sys / 3600)
        m = floor(Int64, (dt_sys % 3600) / 60)
        s = floor(Int64, dt_sys % 60)
        runtime[1] = @sprintf("%02d:%02d:%02d", h, m, s)
        display_values[1] = true
    else
        display_values[1] = false
    end

    return nothing
end

"""
    continue_solver(t::Vector{T}, dt::Vector{T}, tf::T,
                    timer::TimeLimit, show_progress::Bool) where T <: AbstractFloat

Checks whether to continue running the ODE solver. The solver stops if the simulation
finishes, the runtime exceeds the limit set by `timer` or the time step is too small.

Required parameters: `t`, `dt`, `tf`, `timer`, `show_progress`
"""
function continue_solver(t::Vector{T}, dt::Vector{T}, tf::T,
                         timer::TimeLimit, show_progress::Bool) where T <: AbstractFloat
    wtime_min = timer.wtime_min
    time_sys = timer.time_sys

    if !isinf(wtime_min)
        if time() > time_sys[1] + 60*wtime_min
            # note: hack seems to work, but need to maintain # blank lines
            if show_progress
                println("\n\n\n\n")
            else
                println("")
            end
            # note: @code_warntype shows logger and err are type-unstable (from @warn)
            #       logger::Union{Nothing, Base.CoreLogging.AbstractLogger} (red)
            #       err::Any (red)
            # doesn't seem to hurt performance, would rather warn than print
            @warn "Exceeded time limit of $wtime_min minutes (stopping evolve_ode!...)\n"
            return false
        end
    end
    if t[1] + dt[1] == t[1]
        @warn "Time step dt = $(dt[1]) is too small (stopping evolve_ode!...)"
        return false
    end
    return t[1] < tf
end
