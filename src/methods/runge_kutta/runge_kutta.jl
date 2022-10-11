
struct RungeKutta <: DiffEqMethod 
    name::Symbol
    butcher::Matrix{<:AbstractFloat}
    iteration::Iteration
end