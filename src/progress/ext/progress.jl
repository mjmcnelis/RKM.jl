# TODO: make docstring
function RKM.create_progress(n; showspeed = true, color = :gray)
    return Progress(n; showspeed, color)
end

"""
    RKM.monitor_progress(t::Vector{T}, progress::Progress,
                     checkpoints::Vector{T}) where T <: AbstractFloat

Updates the `progress` meter percentage points depending on how many
`checkpoints` the current time `t` has passed.

Required parameters: `t`, `progress`, `checkpoints`
"""
function RKM.monitor_progress(t::Vector{T}, progress::Progress,
                              checkpoints::Vector{T}) where T <: AbstractFloat
    if length(checkpoints) > 1
        dt = checkpoints[2] - checkpoints[1]
        idx = Int(floor(Float64((t[1] - checkpoints[1])/dt))) + 1
        idx = min(length(checkpoints), idx)
        for i in 1:idx
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
