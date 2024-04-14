"""
$(TYPEDEF)

Stores the Butcher tableau and properties of a given Runge-Kutta method.

# Fields
$(TYPEDFIELDS)
"""
struct RungeKutta{T, S, S2, P, I, F} <: ODEMethod where {T <: AbstractFloat, S, S2, P,
                                                         I <: Iteration, F <: Function}
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
    order::SVector{P, T}
    """Determines whether the method is explicit or full/diagonal implicit"""
    iteration::I
    """Determines whether the method has the FSAL property"""
    fsal::Bool
    """Determines which stages are explicit"""
    explicit_stage::SVector{S, Bool}
    """Abbreviated name for the Runge-Kutta method"""
    code_name::String
    """Function used to reconstruct Runge-Kutta method"""
    reconstructor::F
end

# TODO: implement option to use Euler or generic 2nd order
#       should I just replace the embedded row and rename label?

"""
    RungeKutta(name::Symbol, butcher::SMatrix{N, M, T, NM}, iteration::Iteration,
               reconstructor::Function) where {N, M, T <: AbstractFloat, NM}

Outer constructor for `RungeKutta`.

Required parameters: `name`, `butcher`, `iteration`, `reconstructor`
"""
function RungeKutta(name::Symbol, butcher::SMatrix{N, M, T, NM}, iteration::Iteration,
                    reconstructor::Function) where {N, M, T <: AbstractFloat, NM}

    nrow, ncol = size(butcher)
    stages = ncol - 1                               # number of stages
    p = nrow - ncol + 1                             # number of primary/embedded updates

    # split butcher tableau
    c   = butcher[1:ncol-1, 1] |> SVector{stages, T}
    A_T = butcher[1:ncol-1, 2:ncol] |> transpose |> SMatrix{stages, stages, T, stages^2}
    b   = butcher[ncol, 2:ncol] |> SVector{stages, T}
    # TODO: generalize b_hat to multiple embedded pairs
    b_hat = butcher[nrow, 2:ncol] |> SVector{stages, T}

    # properties
    order          = order_prop(name, T, p)         # order of each RK update
    fsal           = method_is_fsal(butcher)        # whether method has FSAL property
    explicit_stage = explicit_stage_prop(butcher)   # mark which stages are explicit

    # code name label
    code_name = make_code_name(name)

    return RungeKutta(name, c, A_T, b, b_hat, stages, order, iteration,
                      fsal, explicit_stage, code_name, reconstructor)
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

function list_explicit_runge_kutta_methods()
    println("Standard (compatible with adaptive = Fixed(), Doubling()):\n")
    """
    \tLow order (1-3)    | Euler1, Heun2, Midpoint2, Ralston2, Generic2, Heun3,
    \t                   | Ralston3, RungeKutta3, ShuOsher3, SpiteriRuuth3, Generic3
    \t-------------------------------------------------------------------------------
    \tMedium order (4-6) | RungeKutta4, ThreeEightsRule4, Ralston4, Ketcheson4,
    \t                   | Butcher5, Butcher6
    \t-------------------------------------------------------------------------------
    \tHigh order (7-9)   | Curtis8, Shanks8, ShanksPseudo8
    """ |> println

    println("\nEmbedded (compatible with adaptive = Fixed(), Embedded(), Doubling()):\n")
    """
    \tLow order (1-3)       | Fehlberg12, HeunEuler21, BogackiShampine32
    \t-------------------------------------------------------------------------------
    \tMedium order (4-6)    | Fehlberg45, CashKarp54, DormandPrince54,
    \t                      | BogackiShampine54, Tsitouras54, Verner56, Verner65
    \t-------------------------------------------------------------------------------
    \tHigh order (7-9)      | Fehlberg78, DormandPrince87
    \t-------------------------------------------------------------------------------
    \tVery high order (10+) | Feagin108, Feagin1210, Feagin1412
    """ |> println
end

function list_implicit_runge_kutta_methods()
    println("Standard (compatible with adaptive = Fixed(), Doubling()):\n")
    """
    \tLow order (1-3)    | BackwardEuler1, ImplicitMidpoint2, TrapezoidRuleBDF2,
    \t                   | QinZhang2, KraaijevangerSpijker2, PareschiRusso2,
    \t                   | PareschiRusso3, Crouzeix3, DIRKL3
    \t-------------------------------------------------------------------------------
    \tMedium order (4-6) | Norsett4
    """ |> println

    println("\nEmbedded (compatible with adaptive = Fixed(), Embedded(), Doubling()):\n")
    """
    \tLow order (1-3)    | CrankNicolson21, LobattoIIIB21
    \t-------------------------------------------------------------------------------
    \tMedium order (4-6) | LobattoIIICS42
    """ |> println
end
