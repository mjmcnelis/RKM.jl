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
    wtime_min::Int64
    """Determines whether or not the runtime has exceeded the time limit"""
    past_time::MVector{10,Bool}
    """Counter for the number of time steps done so far by the solver"""
    counter::MVector{1,Int64}
    # lk
    tsk::Task
end

"""
    function TimeLimit(; wtime_min::Int64 = 60)

Outer constructor for `TimeLimit`.
"""
function TimeLimit(; wtime_min::Int64 = 60)
    past_time = MVector{10,Bool}(repeat([false],10)...)
    counter = MVector{1,Int64}(0)

    timer() = sleep(60 * wtime_min)
    tsk = Task(timer)
    @show  wtime_min
    schedule(tsk)

    return TimeLimit(wtime_min, past_time, counter, tsk)
end

"""
    reset_timer!(timer::TimeLimit)

Resets the `timer` fields to values given by `TimeLimit` outer constructor

Required parameters: `timer`
"""
function reset_timer!(timer::TimeLimit)
    timer.past_time .= false
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
    @unpack wtime_min, past_time, tsk = timer

    # schedule(tsk)
    # yieldto(tsk)
    # schedule(tsk)
    # @show istaskdone(tsk)
    # q()
    # @show past_time

    # Threads.@spawn begin
    return nothing
        #=
    @async begin
        # @info "Timer set to $wtime_min minute(s)"
        sleep(60 * wtime_min)
        # sleep(0)

        while !all(past_time)#islocked(lk)
            past_time .= true
            println(past_time)
        end

        # println("aji")
        past_time[1] = true
    end

        # unlock(task)
        # while !past_time[1]
        #     println(past_time[1])
        #     past_time[1] = true
        # end
        #=
        function a(past_time)
            past_time[1] = true
        end
        t = @task a(past_time)
        # @show "hi"
        # yieldto(t)
        schedule(t);
        # yieldto(t)
        # past_time[1] = true
        =#
    # can I make use of lock/unlock to schedule task
    =#
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
    @unpack counter, wtime_min, past_time, tsk = timer
    # TODO: add counter somewhere else?
    counter[1] += 1
    # if istaskdone(tsk)
    # if past_time[rand(1:10)]
    #     @warn "Exceeded time limit of $wtime_min minutes (stopping evolve_ode!...)\n"
    #     return false
    # end
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
