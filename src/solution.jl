
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
