
# TODO: make docstrings

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