abstract type Magnus <: ODEMethod end

"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
struct MagnusAugmented <: Magnus
    """Name of the Magnus method"""
    name::Symbol
    """Number of stages in the Magnus method"""
    stages::Int64
    """Order of the Magnus update"""
    order::Int64
    """Nonlinear Magnus method is always implicit"""
    iteration::SingleImplicit
    """Abbreviated name for the Magnus method"""
    code_name::String
end

"""
    MagnusAugmented(; order::Int64)

Outer constructor for `MagnusAugmented`.

Required parameters: `order`

"""
function MagnusAugmented(; order::Int64)
    @assert order ∈ [2, 4] "order = $order must be 2 or 4"

    name = Symbol("Magnus_Augmented_$(order)")

    stages = 1

    # note: no need for precision?
    # i need dual jacobian caches, but maybe not here?
    # TODO: is there a more convenient way to group terms? (via matrix)

    iteration = SingleImplicit()

    # code name label
    code_name = make_code_name(name)

    return MagnusAugmented(name, stages, order, iteration, code_name)
end

"""
$(TYPEDEF)

# Fields
$(TYPEDFIELDS)
"""
struct MagnusGenerator <: Magnus
    """Name of the Magnus method"""
    name::Symbol
    """Number of stages in the Magnus method"""
    stages::Int64
    """Order of the Magnus update"""
    order::Int64
     """Nonlinear Magnus method is always implicit"""
    iteration::SingleImplicit
    """Whether to make the Magnus update time symmetric (order increases by 1 if odd)"""
    time_symmetric::Bool
    """Abbreviated name for the Magnus method"""
    code_name::String
end

"""
    MagnusGenerator(; order::Int64, time_symmetric::Bool = false)

Outer constructor for `MagnusGenerator`.

Required parameters: `order`, `time_symmetric`

"""
function MagnusGenerator(; order::Int64, time_symmetric::Bool = false)
    @assert 1 <= order <= 3 "order = $order is not valid"

    name = Symbol("Magnus_Generator_$(order)")

    stages = 1

    # note: no need for precision?
    # i need dual jacobian caches, but maybe not here?
    # TODO: is there a more convenient way to group terms? (via matrix)

    iteration = SingleImplicit()

    # code name label
    code_name = make_code_name(name)

    return MagnusGenerator(name, stages, order, iteration, time_symmetric, code_name)
end

function Base.show(io::IO, M::Magnus)
    for field in M |> typeof |> fieldnames
        println("$field = $(getproperty(M, field))")
    end
end
