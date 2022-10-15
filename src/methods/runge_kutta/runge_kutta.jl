
abstract type RungeKutta <: OrdinaryDiffEqMethod end
struct FixedRungeKutta <: RungeKutta
    # TODO: transpose butcher table? 
    # note: transpose operation is ' (e.g. butcher' .|> precision)
    name::Symbol
    butcher::Matrix{<:AbstractFloat}
    iteration::Iteration
end

struct EmbeddedRungeKutta <: RungeKutta 
    name::Symbol
    butcher::Matrix{<:AbstractFloat}
    iteration::Iteration
    # TODO: implement option to use Euler or generic 2nd order 
    #       should I just replace the embedded row and rename label? 
    # TODO: add field for FSAL
end

function Base.show(io::IO, RK::RungeKutta)
    print(io, "$(str_name(RK.name))")#\n$(DataFrame(RK.butcher, :auto))")
end