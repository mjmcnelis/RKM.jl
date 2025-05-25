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

    set_runtime_display!(timer)
    @unpack runtime, display_values = timer

    generate_showvalues(runtime, t) = () -> [("runtime", runtime[1]),
                                              ("t", t[1]),]

    if length(checkpoints) > 1
        dt_check = checkpoints[2] - checkpoints[1]
        idx = floor(Int64, Float64((t[1] - checkpoints[1])/dt_check)) + 1
        idx = min(length(checkpoints), idx)
        for i in 1:idx
            popfirst!(checkpoints)
        end
        if display_values[1]
            showvalues = generate_showvalues(runtime, t)
            update!(progress, 100 - length(checkpoints); showvalues)
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
