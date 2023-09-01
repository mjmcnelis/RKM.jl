
# this is a minor modification to the ^ operator in DoubleFloats.jl
# to prevent allocations when isinteger(n) == true
# for cases where beta1, beta2, beta3, alpha2 or alpha3 is zero
# original code: https://github.com/JuliaMath/DoubleFloats.jl/blob/87ba8de8f3b6c412f899f02b23f7b9bb9b70f65d/src/math/elementary/explog.jl#L123C1-L123C9
function Base.:(^)(r::DoubleFloat{T}, n::DoubleFloat{T}) where {T <: IEEEFloat}
    if isinteger(n)
        # convert n to Float64 first
        return r^Int(Float64(n))
    else
       return exp(n * log(r))
    end
end