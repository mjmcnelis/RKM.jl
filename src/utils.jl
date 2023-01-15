
# TODO: make docstrings 

function rationalize(x::Float64; sigdigits = 16)
    # TODO: generalize sigdigits for any precision 
    fraction = Int(round(x*10^(sigdigits-1),digits=0))//10^(sigdigits-1)
    return fraction
end

function calloc_vector(static_array::Bool, precision::Type{T}, 
                       rows::Int64) where T <: AbstractFloat
    if static_array
        return @MVector zeros(precision, rows)
    end
    return zeros(precision, rows)
end 

function calloc_matrix(static_array::Bool, precision::Type{T}, 
                       rows::Int64, cols::Int64) where T <: AbstractFloat
    if static_array
        return @MMatrix zeros(precision, rows, cols)
    end
    return zeros(precision, rows, cols)
end 