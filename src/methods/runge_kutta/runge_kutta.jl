
struct RungeKutta <: ODEMethod
    # TODO: transpose butcher table? 
    # note: transpose operation is ' (e.g. butcher' .|> precision)
    name::Symbol
    butcher::Matrix{<:AbstractFloat}
    iteration::Iteration
    fsal::FirstSameAsLast
end

# TODO: implement option to use Euler or generic 2nd order 
    #       should I just replace the embedded row and rename label? 
    # TODO: add field for FSAL

# Properties = Dict(:iteration => Explicit(),
#                   :fsal => FSAL(),
#                   :symp
#                   :)

function RungeKutta(; name::Symbol, butcher::Matrix{<:AbstractFloat})
    # determine properties 
    iteration = iteration_prop(butcher)
    fsal      = fsal_prop(butcher)

    # @show name fsal 
    # println("")

    RungeKutta(name, butcher, iteration, fsal)
end

function Base.show(io::IO, RK::RungeKutta)
    print(io, "$(str_name(RK.name))")#\n$(DataFrame(RK.butcher, :auto))")
end