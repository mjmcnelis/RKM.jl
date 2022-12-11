
struct RungeKutta{T} <: ODEMethod where T <: AbstractFloat
    # TODO: transpose butcher table? 
    # note: transpose operation is ' (e.g. butcher' .|> precision)
    name::Symbol
    c::Vector{T}
    A::Matrix{T}
    b::Vector{T}
    b_hat::Vector{T}
    stages::Int64
    nrow::Int64     # TEMP
    ncol::Int64
    precision::Type{T}
    order::Vector{T}
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
    precision = precision_prop(butcher)         # determine properties
    order     = order_prop(name, butcher)
    iteration = iteration_prop(butcher)
    fsal      = fsal_prop(butcher)
    code_name = make_code_name(name)            # get code name label 

    nrow, ncol = size(butcher)  
    stages = ncol - 1

    c = butcher[1:ncol-1, 1]
    A = butcher[1:ncol-1, 2:ncol]
    b = butcher[ncol, 2:ncol]
    # TODO: change nrow -> ncol + i (where i is the ith embedded pair)
    #       would be necessary when have multiple embedded pairs
    b_hat = butcher[nrow, 2:ncol]
   
    RungeKutta(name, c, A, b, b_hat, stages, nrow, ncol, precision, 
               order, iteration, fsal, code_name)
end

function Base.show(io::IO, RK::RungeKutta)
    print(io, "$(str_name(RK.name))")#\n$(DataFrame(RK.butcher, :auto))")
end