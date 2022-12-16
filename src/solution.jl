
# TODO: we could also try sol.y = DataFrame
#       not sure if it's difficult to append a vector (not DataFrame)

abstract type DataFormat end 
struct TimeSlice <: DataFormat end      # better for large systems
struct SpaceSlice <: DataFormat end       # better for small systems

# y(0) = [1, 2]
# y(1) = [3, 4]
# push (Time)
# sol.y = [[1,2], [3,4]]
# append (Space) 
# sol.y = [[1,3], [2,4]]

struct Solution{T1 <: Vector{<:Vector{<:AbstractFloat}}, T2 <: Vector{<:AbstractFloat},
                T3 <: DataFormat}
    y::T1 
    t::T2
    data_format::T3
end

function Solution(; precision::Type{<:AbstractFloat} = Float64, dimensions::Int64,
                    # TODO: remove default value 
                    data_format::DataFormat = SpaceSlice()) 

    y = data_format isa SpaceSlice ? Vector{Vector{precision}}([[] for i = 1:dimensions]) :
                                     Vector{Vector{precision}}()
    
    t = Vector{precision}()
    Solution(y, t, data_format)
end

function update_solution!(sol, y::Vector{T}, t::MVector{1,T}) where T <: AbstractFloat
    # TODO: do type dispatching
    if sol.data_format isa TimeSlice 
        push!(sol.y, copy(y))
    elseif sol.data_format isa SpaceSlice 
        for i in eachindex(y) 
            append!(sol.y[i], y[i])
        end
    end
    append!(sol.t, t[1])
    nothing 
end

function sizehint_solution!(sol, t_span)
    @unpack t0, tf, dt0 = t_span
    steps = round((tf - t0)/dt0) |> Int64

    # TODO: what to do for TimeSlice format? 

    # note: this is only useful if we're saving at regular intervals
    for i in eachindex(sol.y)
        # TODO: why is steps + 1 not sufficient?
        sizehint!(sol.y[i], steps + 2)
    end
    sizehint!(sol.t, steps + 2)
    nothing
end
