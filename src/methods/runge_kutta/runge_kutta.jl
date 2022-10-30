
struct RungeKutta{T1} <: ODEMethod
    # TODO: transpose butcher table? 
    # note: transpose operation is ' (e.g. butcher' .|> precision)
    name::Symbol
    butcher::Matrix{T1}
    # TODO: should I just wrap this in a Properties struct? 
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

    RungeKutta(name, butcher, iteration, fsal)
end

function Base.show(io::IO, RK::RungeKutta)
    print(io, "$(str_name(RK.name))")#\n$(DataFrame(RK.butcher, :auto))")
end