
struct Solution{T1 <: Vector{<:Vector{<:AbstractFloat}}, T2 <: Vector{<:AbstractFloat}}
    y::T1 
    t::T2
end

function Solution(; precision::Type{<:AbstractFloat} = Float64)
    y = Vector{Vector{precision}}() 
    t = Vector{precision}()
    Solution(y, t)
end

