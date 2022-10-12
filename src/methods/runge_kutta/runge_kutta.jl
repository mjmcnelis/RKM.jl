
struct RungeKutta <: DiffEqMethod 
    # TODO: transpose butcher table? 
    # note: transpose operation is ' (e.g. butcher' .|> precision)
    name::Symbol
    butcher::Matrix{<:AbstractFloat}
    iteration::Iteration
    # TODO: add Fixed(), Embedded()?
    # Doubling, FiniteDiff can only use Fixed RK schemes
    # Embedded can only use Embedded RK schemes
end

function Base.show(io::IO, RK::RungeKutta)
    print(io, "$(str_name(RK.name))")#\n$(DataFrame(RK.butcher, :auto))")
end