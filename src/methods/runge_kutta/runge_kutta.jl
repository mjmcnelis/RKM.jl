
struct RungeKutta <: DiffEqMethod 
    name::Symbol
    butcher::Matrix{<:AbstractFloat}
    iteration::Iteration
end

function Base.show(io::IO, RK::RungeKutta)
    println(io, "$(str_name(RK.name))\n$(DataFrame(RK.butcher, :auto))")
end