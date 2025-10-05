"""
$(TYPEDEF)

Stores the Butcher tableau and properties of a given Runge-Kutta method.

# Fields
$(TYPEDFIELDS)
"""
struct RungeKutta{T, S, S2, E, SE,
                  R, Q, RQ, P, I, F} <: ODEMethod where {T <: AbstractFloat,
                                                         S, S2, E, SE, R, Q, RQ, P,
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
    b_hat::SMatrix{E, S, T, SE}
    """Polynomial coefficients used for continuous output"""
    ω::SMatrix{R, Q, T, RQ}
    """Number of stages in the Runge-Kutta method"""
    stages::Int64
    """Order of the primary (and embedded) update(s) of the Runge-Kutta method"""
    order::SVector{P, T}
    """Determines whether the method is explicit or full/diagonal implicit"""
    iteration::I
    """Determines whether the method has the FESAL property (see ../properties/fesal.jl)"""
    fesal::Bool
    """Determines which stages are explicit"""
    explicit_stage::SVector{S, Bool}
    """Abbreviated name for the Runge-Kutta method"""
    code_name::String
    """Function used to reconstruct Runge-Kutta method"""
    reconstructor::F
end

"""
    RungeKutta(name::Symbol, butcher::SMatrix{N, M, T, NM},
               iteration::Iteration, reconstructor::Function;
               ω::SMatrix{R, Q, T, RQ} = SMatrix{0,0,T,0}()
              ) where {T <: AbstractFloat, N, M, NM, R, Q, RQ}

Outer constructor for `RungeKutta`.

Required parameters: `name`, `butcher`, `iteration`, `reconstructor`

Note: `ω` are interpolation coefficients used for continuous output,
      but some methods may not have them.
"""
function RungeKutta(name::Symbol, butcher::SMatrix{N, M, T, NM},
                    iteration::Iteration, reconstructor::Function;
                    ω::SMatrix{R, Q, T, RQ} = SMatrix{0,0,T,0}()
                   ) where {T <: AbstractFloat, N, M, NM, R, Q, RQ}

    nrow, ncol = size(butcher)
    stages = ncol - 1                               # number of stages
    p = nrow - ncol + 1                             # number of primary/embedded updates

    # split butcher tableau
    c = butcher[1:ncol-1, 1] |> SVector{stages, T}
    A_T = butcher[1:ncol-1, 2:ncol] |> transpose |> SMatrix{stages, stages, T, stages^2}
    b = butcher[ncol, 2:ncol] |> SVector{stages, T}
    b_hat = butcher[ncol+1:end, 2:ncol] |> SMatrix{p-1, stages, T, (p-1)*stages}

    # TODO: do I really need to pipeline it like this? if I don't,
    #       then ω appears as Core.Const(...) in @code_warntype
    ω = Matrix(ω) |> SMatrix{R, Q, T, RQ}

    # get properties
    order = order_prop(name, T, p)                  # order of each RK update
    fesal = get_fesal(A_T, b, c)                    # whether method has FESAL property
    explicit_stage = explicit_stage_prop(butcher)   # mark which stages are explicit

    # code name label
    code_name = make_code_name(name)

    return RungeKutta(name, c, A_T, b, b_hat, ω, stages, order, iteration,
                      fesal, explicit_stage, code_name, reconstructor)
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
    """
    Low order (1-3)       | Euler1, Heun2, Midpoint2, Ralston2, Fehlberg2, Heun2, Heun3,
                          | Ralston3, Kutta3, ShuOsher3, SpiteriRuuth3, BogackiShampine3
    -------------------------------------------------------------------------------
    Medium order (4-6)    | RungeKutta4, ThreeEightsRule4, Ralston4, Ketcheson4, Butcher5,
                          | Fehlberg5, CashKarp5, DormandPrince5, BogackiShampine5,
                          | Tsitouras5, Verner5, Butcher6, Verner6
    -------------------------------------------------------------------------------
    High order (7-9)      | Fehlberg7, DormandPrince8, Curtis8, Shanks8, ShanksPseudo8
    -------------------------------------------------------------------------------
    Very high order (10+) | Feagin10, Feagin12, Feagin14
    """ |> println

    return nothing
end

function list_diagonal_implicit_runge_kutta_methods()
    """
    Low order (1-3)       | BackwardEuler1, TrapezoidRuleBDF2, ImplicitTrapezoid2,
                          | ImplicitMidpoint2, QinZhang2, KraaijevangerSpijker2,
                          | PareschiRusso2, LobattoIIIB2, PareschiRusso3, Crouzeix3,
                          | DIRKL3
    -------------------------------------------------------------------------------
    Medium order (4-6)    | Norsett4, LobattoIIICS42
    """ |> println

    return nothing
end
