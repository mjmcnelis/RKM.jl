
struct LinearMultistep{T, S, P} <: ODEMethod where {T <: AbstractFloat, S, P}
    name::Symbol
    b::SVector{S,T}
    stages::Int64
    order::SVector{P,T}
    code_name::String
end

function LinearMultistep(; name::Symbol, table::Vector{T}) where T <: AbstractFloat
    order     = order_prop(name, T)
    code_name = make_code_name(name)

    ncol, = size(table)
    stages = ncol - 1

    b = table[2:ncol] |> SVector{stages}

    return LinearMultistep(name, b, stages, order, code_name)
end

function Base.show(io::IO, LM::LinearMultistep)
    for field in LM |> typeof |> fieldnames
        println("$field = $(getproperty(LM, field))")
    end
end