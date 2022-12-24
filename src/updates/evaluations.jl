
function add_function_evaluations!(::Explicit, ::Fixed, FE::MVector{1,Int64},
                                   method::RungeKutta, args...)
    @.. FE += method.stages
    nothing
end

function add_function_evaluations!(::Explicit, ::Doubling, FE::MVector{1,Int64},
                                   method::RungeKutta, attempts::Int64)
    evals = 1 + attempts*(3*method.stages - 2)
    @.. FE += evals
    nothing
end

function add_function_evaluations!(::Explicit, ::Embedded, FE::MVector{1,Int64},
                                   method::RungeKutta, attempts::Int64)
    # TODO: have not implemented FSAL yet
    evals = 1 + attempts*(method.stages - 1)
    @.. FE += evals
    nothing
end