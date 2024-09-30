
create_progress(args...; kwargs...) = nothing

function monitor_progress(args...)
    @warn "Progress meter will not display unless RKMProgressMeterExt is activated. \
           Install ProgressMeter in your base environment and load the module." maxlog=1
    return nothing
end