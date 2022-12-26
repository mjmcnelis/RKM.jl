
# y(0) = [1, 2]
# y(1) = [3, 4]
# y(2) = [5, 6]
# append! (LinearColumn)
# sol.y = [1,2,3,4,5,6]

struct Solution{T <: AbstractFloat}
    y::Vector{T}
    t::Vector{T}
    # TODO: separate from y(t) solution
    FE::MVector{1,Int64}
    dimensions::Int64
end

function Solution(; precision::Type{<:AbstractFloat}, dimensions::Int64)
    y = Vector{precision}()                                 
    t = Vector{precision}()
    FE = MVector{1,Int64}(0)

    Solution(y, t, FE, dimensions)
end

function update_solution!(sol, y::Vector{T}, t::MVector{1,T2}) where {T <: AbstractFloat, 
                                                                      T2 <: AbstractFloat}
    append!(sol.y, y)
    append!(sol.t, t[1])
    nothing 
end

# TODO: whether or not I call this depends if save at regular intervals
function sizehint_solution!(sol, t_span, dimensions)
    @unpack t0, tf, dt0 = t_span
    steps = round((tf - t0)/dt0) |> Int64

    sizehint!(sol.y, dimensions*(steps + 2))
    sizehint!(sol.t, steps + 2)
    nothing
end

function reshape_solution(sol)
    y, t = sol.y, sol.t
    # TODO: replace length(t) if use deleteat for PDEs
    y = reshape(y, sol.dimensions, length(t))'
    y, t
end
