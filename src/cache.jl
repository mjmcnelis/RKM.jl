
abstract type RKMCache end

struct UpdateCache{T <: AbstractFloat} <: RKMCache
    dy::Matrix{T}
    dy_LM::Matrix{T}
    y::Vector{T}
    y_tmp::Vector{T}
    f_tmp::Vector{T}
    f::Vector{T}
    J::Matrix{T}
    y1::Vector{T}
    y2::Vector{T}
    error::Vector{T}
end

function UpdateCache(precision::Type{T}, y::Vector{T}, method::ODEMethod,
                     adaptive::AdaptiveStepSize,
                     dimensions::Int64) where T <: AbstractFloat

    @unpack iteration, stages = method

    n = iteration isa Explicit ? 0 : dimensions
    m = adaptive isa Fixed ? 0 : dimensions
    p = iteration isa Explicit && adaptive isa Fixed ? 0 : dimensions

    if method isa LinearMultistep
        @unpack start_method = method
        dy = zeros(precision, dimensions, start_method.stages)
        dy_LM = zeros(precision, dimensions, stages)
    else
        dy = zeros(precision, dimensions, stages)
        dy_LM = Array{precision}(undef, 0, 0)
    end

    y_tmp = zeros(precision, dimensions)
    f_tmp = zeros(precision, dimensions)
    f     = zeros(precision, dimensions)
    # TODO: should have option to make it sparse
    J     = zeros(precision, n, n)
    y1    = zeros(precision, m)
    y2    = zeros(precision, m)
    error = zeros(precision, p)

    return UpdateCache(dy, dy_LM, y, y_tmp, f_tmp, f, J, y1, y2, error)
end
