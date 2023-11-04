
function solve_linear_tmp!(cache::LinearCache)
    # TODO: extend to methods that don't involve factorization
    if cache.isfresh
        # note: allocates on do_factorization (any way to avoid this?)
        fact = do_factorization_tmp(cache, cache.alg)
        # note: cache.isfresh is set to false by set_cacheval
        cache.cacheval = fact
        cache.isfresh = false
    end
    _ldiv!(cache.u, cache.cacheval, cache.b)
    return nothing
end

function do_factorization_tmp(cache, alg::AbstractFactorization)
    do_factorization(alg, cache.A, cache.b, cache.u)
end

# note: much simpler than https://github.com/SciML/LinearSolve.jl/blob/main/src/factorization.jl
# ln.353-372 so be wary
_do_factorization(cache, A::SparseMatrixCSC) = lu!(cache.cacheval, A)
function do_factorization_tmp(cache, ::KLUFactorization)
    # note: really UMFPACKFactorization
    _do_factorization(cache, cache.A)
end


