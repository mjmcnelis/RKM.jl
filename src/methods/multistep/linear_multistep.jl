
abstract type LinearMultistep <: ODEMethod end

struct Adams{T, S} <: LinearMultistep where {T <: AbstractFloat, S}
    name::Symbol
    b::SVector{S,T}
    b_pred::SVector{S,T}
    stages::Int64
    order::Int64
    iteration::Iteration
    code_name::String
end

function Adams(; name::Symbol, order::Int64, table::Matrix{T},
                 table_pred::Matrix{T} = Array{Float64}(undef, 0, 0)
              ) where T <: AbstractFloat

    code_name = make_code_name(name)
    stages = order

    b = table[order, 2:order+1] |> SVector{order}

    if isempty(table_pred)
        b_pred = b  # dummy filler
        iteration = Explicit()
    else
        b_pred = table_pred[order, 2:order+1] |> SVector{order}
        iteration = SingleImplicit()
    end

    # TODO: replace stages with steps
    return Adams(name, b, b_pred, stages, order, iteration, code_name)
end

function Base.show(io::IO, LM::LinearMultistep)
    for field in LM |> typeof |> fieldnames
        println("$field = $(getproperty(LM, field))")
    end
end