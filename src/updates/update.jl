
@inline function sensitivity_update!(sensitivity::NoSensitivity,
                                     update_cache::UpdateCache, method::ODEMethod)
    return nothing
end

@inline function sensitivity_update!(sensitivity::DecoupledDirect,
                                     update_cache::UpdateCache, method::ODEMethod)

    S = update_cache.S
    S_tmp = update_cache.S_tmp
    dS = update_cache.dS

    b = method.b
    stages = method.stages

    @.. S_tmp = S
    for j in 1:stages
        if iszero(b[j])
            continue
        end
        @.. S_tmp = S_tmp + b[j]*dS[:,:,j]
    end

    return nothing
end