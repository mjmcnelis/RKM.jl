
create_progress(args...; kwargs...) = nothing

function monitor_progress(args...)
    @warn "Progress meter will not display unless RKMProgressMeterExt is activated. \
           Install and load the module ProgressMeter in your base environment." maxlog=1
    return nothing
end