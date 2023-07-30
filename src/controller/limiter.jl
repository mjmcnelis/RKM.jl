
abstract type LimiterMethod end

# TODO: make piecewise and smooth limiter methods
# eventually these would also impose bounds by dt_min, dt_max
# if I re-evaluate rescale based on dt_min <= dt <= dt_max,
# would floating precision errors be a problem?

function limit_time_step(controller::PIDControl, rescale)
    @unpack safety, low, high = controller
    return min(high, max(low, safety*rescale))
end