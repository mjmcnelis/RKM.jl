
abstract type RootMethod end 
struct Newton <: RootMethod end         # TODO: add line search later
struct FixedPoint <: RootMethod end

abstract type StageFinder end 

# good enough start (wrap caches later)
# TODO: determine Jacobian method here (finite or forward diff)
@kwdef struct ImplicitStageFinder <: StageFinder
    root_method::RootMethod = Newton()
    epsilon::Float64        = 1.0e-8    # TODO: reuse adaptive epsilon or 100x smaller?
    max_iterations::Int64   = 10
end
