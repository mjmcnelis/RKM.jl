
@kwdef struct TimeSpan
    t0::Float64
    tf::Float64
    dt0::Float64
end

struct TimeLimit
    wtime_min::Int64
    time_limit::DateTime
    frequency::Int64
    counter::Vector{Int64}
end

function TimeLimit(; wtime_min = 60, frequency = 100)
    time_limit = Dates.now() + Dates.Minute(round(wtime_min))
    counter = [0]

    TimeLimit(wtime_min, time_limit, frequency, counter)
end

function check_time(t::MVector{1,Float64}, tf::Float64, timer::TimeLimit)
    t[1] < tf && !past_time_limit(timer)
end

function past_time_limit(timer)
    @unpack wtime_min, counter, frequency, time_limit = timer

    if (counter[1] += 1) % frequency == 0 && Dates.now() > time_limit
        @warn "\nExceeded time limit of $(wtime_min) minutes (stopping evolve)\n"
        true
    else
        false
    end
end
