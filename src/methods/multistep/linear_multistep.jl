
struct LinearMultistep{T, S, P} <: ODEMethod where {T <: AbstractFloat, S, P}
    name::Symbol
    b::SVector{S,T}
    b_pred::SVector{S,T}
    stages::Int64
    order::SVector{P,T}
    iteration::Iteration
    code_name::String
end

function LinearMultistep(; name::Symbol, table::Matrix{T}) where T <: AbstractFloat
    order     = order_prop(name, T)
    code_name = make_code_name(name)

    nrow, ncol = size(table)
    @assert 1 <= nrow <= 2 && ncol >= 2

    stages = ncol - 1

    b = table[1, 2:ncol] |> SVector{stages}
    b_pred = table[max(1,nrow), 2:ncol] |> SVector{stages}

    iteration = nrow == 1 ? Explicit() : SingleImplicit()

    return LinearMultistep(name, b, b_pred, stages, order, iteration, code_name)
end

function Base.show(io::IO, LM::LinearMultistep)
    for field in LM |> typeof |> fieldnames
        println("$field = $(getproperty(LM, field))")
    end
end