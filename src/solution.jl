
abstract type DataFormat end 
struct TimeSlice <: DataFormat end      # better for large systems
struct SpaceSlice <: DataFormat end     # better for small systems

# y(0) = [1, 2]
# y(1) = [3, 4]
# y(2) = [5, 6]
# push! (Time)
# sol.y = [[1,2], [3,4], [5,6]]
# append! (Space) 
# sol.y = [[1,3,5], [2,4,6]]

struct Solution{T1 <: Vector{<:Vector{<:AbstractFloat}}, 
                T2 <: Vector{<:AbstractFloat},
                T3 <: DataFormat}
    y::T1 
    t::T2
    data_format::T3
    # TODO: separate from y(t) solution
    FE::MVector{1,Int64}
end

function Solution(; precision::Type{<:AbstractFloat}, dimensions::Int64,
                    data_format::DataFormat)

    y = data_format isa SpaceSlice ? Vector{Vector{precision}}([[] for i = 1:dimensions]) :
                                     Vector{Vector{precision}}()
    t = Vector{precision}()

    FE = MVector{1,Int64}(0)

    Solution(y, t, data_format, FE)
end

update_state!(::TimeSlice, sol, y) = push!(sol.y, copy(y))
function update_state!(::SpaceSlice, sol, y)
    for i in eachindex(y) 
        append!(sol.y[i], y[i])
    end
end

function update_solution!(sol, y::Vector{T}, t::MVector{1,T2}) where {T <: AbstractFloat, 
                                                                      T2 <: AbstractFloat}
    update_state!(sol.data_format, sol, y)
    append!(sol.t, t[1])
    nothing 
end

sizehint_state!(::TimeSlice, sol, N) = sizehint!(sol.y, N + 2)
function sizehint_state!(::SpaceSlice, sol, N)
    for i in eachindex(sol.y)
        # TODO: why is steps + 1 not sufficient?
        sizehint!(sol.y[i], N + 2)
    end
end

# TODO: whether or not I call this depends if save at regular intervals
function sizehint_solution!(sol, t_span)
    @unpack t0, tf, dt0 = t_span
    steps = round((tf - t0)/dt0) |> Int64

    sizehint_state!(sol.data_format, sol, steps)
    sizehint!(sol.t, steps + 2)
    nothing
end
