
"""
    monitor_progress(progress::Progress, checkpoints::Vector{T},
                     t::Vector{T}, timer::TimeLimit, dt::Vector{T})

Updates the `progress` bar percentage points if the current time `t` has passed any
`checkpoints`. ODE variables and solver stats are also displayed in real time.

Required parameters: `progress`, `checkpoints`, `t`, `timer`, `dt`
"""
function monitor_progress(progress::Progress, checkpoints::Vector{T}, t::Vector{T},
                          timer::TimeLimit, dt::Vector{T}) where T <: AbstractFloat

    set_runtime_display!(timer)

    runtime = timer.runtime
    display_values = timer.display_values
    total_steps = timer.total_steps

    generate_showvalues(runtime, total_steps, t, dt) = () -> [
        (:runtime, runtime[1]),
        (:total_steps, total_steps[1]),
        (:t, t[1]),
        (:dt, dt[1]),]

    if length(checkpoints) > 1
        dt_check = checkpoints[2] - checkpoints[1]
        idx = floor(Int64, Float64((t[1] - checkpoints[1])/dt_check)) + 1
        idx = min(length(checkpoints), idx)
        for i in 1:idx
            popfirst!(checkpoints)
        end
        if display_values[1]
            showvalues = generate_showvalues(runtime, total_steps, t, dt)
            update!(progress, 100 - length(checkpoints); showvalues)
        end
    else
        if t[1] >= checkpoints[1]
            showvalues = generate_showvalues(runtime, total_steps, t, dt)
            update!(progress, 100; showvalues)
            sleep(1e-3)
        end
    end
    return nothing
end
