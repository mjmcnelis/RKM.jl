
@kwdef struct TimeSpan
    t0::Float64
    tf::Float64
    dt0::Float64
    dt_min::Float64 = 0.0
    dt_max::Float64 = Inf
end
