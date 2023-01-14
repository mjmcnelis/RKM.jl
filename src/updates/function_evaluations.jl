"""
    add_function_evaluations(FE::MVector{1,Int64}, iteration::Iteration, 
                             adaptive::AdaptiveStepSize, ode_method::ODEMethod, 
                             args...)

Add number of function evaluations per time step to `FE`.

Required parameters: `FE`, `iteration`, `adaptive`, `method`
"""
function add_function_evaluations!(FE::MVector{1,Int64}, iteration::Iteration, 
                                   adaptive::AdaptiveStepSize, method::ODEMethod, 
                                   args...)
    @error "No function for $iteration, $adaptive method $(method.name)"
end

function add_function_evaluations!(FE::MVector{1,Int64}, ::Explicit, ::Fixed, 
                                   method::RungeKutta, args...)
    @.. FE += method.stages
    return nothing
end

function add_function_evaluations!(FE::MVector{1,Int64}, ::Explicit, ::Doubling,
                                   method::RungeKutta, attempts::Int64)
    evals = 1 + attempts*(3*method.stages - 2)
    @.. FE += evals
    return nothing
end

function add_function_evaluations!(FE::MVector{1,Int64}, ::Explicit, ::Embedded,
                                   method::RungeKutta, attempts::Int64)
    # TODO: have not implemented FSAL yet
    @unpack stages, fsal = method
    fsal_stage = fsal isa FSAL ? 1 : 0
   
    # note: fsal advantage does not count for rejected attempts
    evals = 1 - fsal_stage + attempts*(stages - 1)
    @.. FE += evals
    return nothing
end
