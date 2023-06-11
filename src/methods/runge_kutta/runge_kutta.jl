"""
$(TYPEDEF)

Stores the Butcher tableau and properties of a given Runge-Kutta method.

# Fields
$(TYPEDFIELDS)
"""
struct RungeKutta{T, S, S2} <: ODEMethod where {T <: AbstractFloat, S, S2}
    """Name of the Runge-Kutta method"""
    name::Symbol
    """Intermediate time update coefficient of each stage in the Butcher tableau"""
    c::SVector{S, T}
    """Transposed intermediate state update coefficients of each stage in the Butcher tableau"""
    A_T::SMatrix{S, S, T, S2}       # TODO: would this be slow for high order methods?
    """Primary state update coefficients of the Butcher tableau"""
    b::SVector{S, T}
    """Embedded state update coefficients (if any) of the Butcher tableau"""
    b_hat::SVector{S, T}            # TODO: generalize b_hat
    """Number of stages in the Runge-Kutta method"""
    stages::Int64
    """Order of the primary (and embedded) update(s) of the Runge-Kutta method"""
    order::Vector{T}
    # TODO: should I just wrap this in a Properties struct?
    """Determines whether the method is explicit or implicit"""
    iteration::Iteration
    """Determines whether the method has the FSAL property"""
    fsal::FirstSameAsLast
    """Determines which stages are explicit"""
    explicit_stage::SVector{S, Bool}
    """Abbreviated name for the Runge-Kutta method"""
    code_name::String
end

# TODO: implement option to use Euler or generic 2nd order
    #       should I just replace the embedded row and rename label?
    # TODO: add field for FSAL

# Properties = Dict(:iteration => Explicit(),
#                   :fsal => FSAL(),
#                   :symp
#                   :)

"""
    RungeKutta(; name::Symbol, butcher::Matrix{T}) where T <: AbstractFloat

Outer constructor for `RungeKutta`. 

Required parameters: `name`, `butcher`
"""
function RungeKutta(; name::Symbol, butcher::Matrix{T}) where T <: AbstractFloat
    order          = order_prop(name, butcher)      # determine properties
    iteration      = iteration_prop(butcher)
    fsal           = fsal_prop(butcher)
    explicit_stage = explicit_stage_prop(butcher)
    code_name      = make_code_name(name)           # get code name label

    nrow, ncol = size(butcher)
    stages = ncol - 1

    # convert butcher tableau into static arrays
    c   = butcher[1:ncol-1, 1] |> SVector{stages}
    A_T = butcher[1:ncol-1, 2:ncol] |> transpose |> SMatrix{stages, stages}
    b   = butcher[ncol, 2:ncol] |> SVector{stages}
    # TODO: generalize b_hat to multiple embedded pairs
    b_hat = butcher[nrow, 2:ncol] |> SVector{stages}

    return RungeKutta(name, c, A_T, b, b_hat, stages, order, 
                      iteration, fsal, explicit_stage, code_name)
end

function Base.show(io::IO, RK::RungeKutta)
    for field in RK |> typeof |> fieldnames
        # display A rather than A_T
        if field == :A_T
            print("A = ")
            display(RK.A_T')
        else
            println("$field = $(getproperty(RK, field))")
        end
    end
end
