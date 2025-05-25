# TODO: make docstring
function RKM.create_progress(n; showspeed = true, color = :gray)
    return Progress(n; showspeed, color)
end

"""
    RKM.monitor_progress(t::Vector{T}, progress::Progress, checkpoints::Vector{T},
                         timer::TimeLimit) where T <: AbstractFloat

Updates the `progress` meter percentage points depending on how many
`checkpoints` the current time `t` has passed.

Required parameters: `t`, `progress`, `checkpoints`, `timer`
"""
function RKM.monitor_progress(t::Vector{T}, progress::Progress, checkpoints::Vector{T},
                              timer::TimeLimit) where T <: AbstractFloat

    @unpack time_sys, runtime, t_prev, total_steps = timer

    dt_sys = floor(Int64, time_sys[2] - time_sys[1])

    if dt_sys > t_prev[1]
        set_runtime!(timer)
        display_values = true
    else
        display_values = false
    end

    generate_showvalues(runtime, t) = () -> [("runtime", runtime[1]),
                                              ("t", t[1]),]

    if length(checkpoints) > 1
        dt_chkpt = checkpoints[2] - checkpoints[1]

        idx = floor(Int64, Float64((t[1] - checkpoints[1])/dt_chkpt)) + 1
        idx = min(length(checkpoints), idx)

        if idx > 0
            @unpack counter = progress
            showvalues = generate_showvalues(runtime, t)
            update!(progress, counter + idx; showvalues)
            for i in 1:idx
                popfirst!(checkpoints)
            end
        elseif display_values
            showvalues = generate_showvalues(runtime, t)
            update!(progress; showvalues)
        end
    else
        if t[1] >= checkpoints[1]
            showvalues = generate_showvalues(runtime, t)
            update!(progress, 100; showvalues)
            sleep(1e-3)
        end
    end
    return nothing
end
