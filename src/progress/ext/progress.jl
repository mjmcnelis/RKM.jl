# TODO: make docstring
function RKM.create_progress(n; showspeed = true, color = :gray)
    return Progress(n; showspeed, color)
end

"""
    RKM.monitor_progress(progress::Progress, checkpoints::Vector{T},
                         t::Vector{T}, timer::TimeLimit, dt::Vector{T})

Updates a progress meter every second in real time. The percentage points
indicate how much time evolution the solver has completed. ODE variables
and solver stats are also displayed in real time.

Required parameters: `progress`, `checkpoints`, `t`, `timer`, `dt`
"""
function RKM.monitor_progress(progress::Progress, checkpoints::Vector{T}, t::Vector{T},
                              timer::TimeLimit, dt::Vector{T}) where T <: AbstractFloat

    set_runtime_display!(timer)
    @unpack runtime, display_values, total_steps = timer

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
