
function explicit_stage_prop(butcher::SMatrix{N, M, T, NM}) where {N, M,
                                                                   T <: AbstractFloat,
                                                                   NM}
    ncol = size(butcher, 2)
    A = butcher[1:(ncol-1), 2:end]
    stages = ncol-1
    explicit_stage = Vector{Bool}(undef, stages)
    for i in 1:stages
        explicit_stage[i] = all(x -> iszero(x), A[i,i:end])
    end
    return explicit_stage |> SVector{stages}
end