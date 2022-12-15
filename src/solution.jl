
struct Solution{T1 <: Vector{<:Vector{<:AbstractFloat}}, T2 <: Vector{<:AbstractFloat}}
    y::T1 
    t::T2
end

function Solution(; precision::Type{<:AbstractFloat} = Float64, dimensions::Int64)
    # y = Vector{Vector{precision}}() 
    y = Vector{Vector{precision}}([[] for i = 1:dimensions]) 
    t = Vector{precision}()
    Solution(y, t)
end

function update_solution!(sol, y::Vector{T}, t::MVector{1,T}) where T <: AbstractFloat
    # push!(sol.y, copy(y))
    for i in eachindex(y) 
        append!(sol.y[i], y[i])
    end
    append!(sol.t, t[1])
    nothing 
end

function sizehint_solution!(sol, t_span)
    @unpack t0, tf, dt0 = t_span
    steps = round((tf - t0)/dt0) |> Int64

    # note: this is only useful if we're saving at regular intervals
    for i in eachindex(sol.y)
        # TODO: why is steps + 1 not sufficient?
        sizehint!(sol.y[i], steps + 2)
    end
    sizehint!(sol.t, steps + 2)
    nothing
end
