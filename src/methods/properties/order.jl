
# note: needed @inline macro to make SVector return type-stable
@inline function order_prop(name::Symbol, precision::Type{T},
                            p::Int) where T <: AbstractFloat
    order = filter.(isdigit, split(string(name), "_"))
    order = parse.(Float64, filter(x -> x != "", order)) .|> precision
    return SVector{p, T}(order)
end