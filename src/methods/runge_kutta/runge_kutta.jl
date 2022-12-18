
struct RungeKutta{T, S, S2} <: ODEMethod where {T <: AbstractFloat, S, S2}
    name::Symbol
    c::SVector{S, T}
    A_T::SMatrix{S, S, T, S2}
    b::SVector{S, T}
    b_hat::SVector{S, T}
    stages::Int64
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

    # convert butcher tableau into static arrays
    c   = butcher[1:ncol-1, 1] |> SVector{stages}
    A_T = butcher[1:ncol-1, 2:ncol] |> transpose |> SMatrix{stages, stages}
    b   = butcher[ncol, 2:ncol] |> SVector{stages}
    # TODO: generalize b_hat to multiple embedded pairs
    b_hat = butcher[nrow, 2:ncol] |> SVector{stages}
   
    RungeKutta(name, c, A_T, b, b_hat, stages, precision, 
               order, iteration, fsal, code_name)
end

function Base.show(io::IO, RK::RungeKutta)
    print(io, "$(str_name(RK.name))")#\n$(DataFrame(RK.butcher, :auto))")
end