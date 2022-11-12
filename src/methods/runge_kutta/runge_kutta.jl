
struct RungeKutta{T1 <: AbstractFloat, 
                  T2 <: AbstractFloat, 
                  T3 <: AbstractFloat} <: ODEMethod
    # TODO: transpose butcher table? 
    # note: transpose operation is ' (e.g. butcher' .|> precision)
    name::Symbol
    butcher::Matrix{T1}
    precision::Type{T2}
    order::Vector{T3}
    # TODO: should I just wrap this in a Properties struct? 
    iteration::Iteration
    fsal::FirstSameAsLast
    code_name::String
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
    precision = precision_prop(butcher)
    order     = order_prop(name, butcher)
    iteration = iteration_prop(butcher)
    fsal      = fsal_prop(butcher)

    # get code name label 
    code_name = make_code_name(name) 
   
    RungeKutta(name, butcher, precision, order, iteration, fsal, code_name)
end

function Base.show(io::IO, RK::RungeKutta)
    print(io, "$(str_name(RK.name))")#\n$(DataFrame(RK.butcher, :auto))")
end