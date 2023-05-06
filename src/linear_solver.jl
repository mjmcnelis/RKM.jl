
function solve_linear_tmp(cache::LinearCache)
    # TODO: extend to methods that don't involve factorization
    if cache.isfresh
        # note: allocates on do_factorization
        fact = do_factorization(cache.alg, cache.A, cache.b, cache.u)
        # note: cache.isfresh is set to false by set_cacheval
        cache = set_cacheval(cache, fact)
    end
    _ldiv!(cache.u, cache.cacheval, cache.b)
    return cache
end
