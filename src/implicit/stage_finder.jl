
abstract type StageFinder end

# good enough start (wrap caches later)
@kwdef struct ImplicitStageFinder{JM, AF, EMM} <: StageFinder where {JM <: JacobianMethod,
                                                               AF <: AbstractFactorization,
                                                               EMM <: EigenMaxMethod}
    root_method::RootMethod = Newton()  # TODO: try RM
    state_jacobian::JM = FiniteJacobian()
    linear_method::AF = LUFactorization()
    eigenmax::EMM = NoEigenMax()
    # TODO: reuse adaptive epsilon or 100x smaller?
    epsilon::Float64 = 1e-8
    max_iterations::Int64 = 10
    # TODO: make outer constructor to check p_norm value
    p_norm::Float64 = 2.0
    # add iterations_per_stage
end

function reconstruct_stage_finder(stage_finder::ImplicitStageFinder,
                                  ode_wrap!::ODEWrapperState, f::Vector{T},
                                  y::Vector{T}) where T <: AbstractFloat
    @unpack state_jacobian, eigenmax = stage_finder

    # TODO: move state_jacobian out of stage_finder?
    state_jacobian = reconstruct_jacobian_method(state_jacobian, ode_wrap!, f, y)
    @set! stage_finder.state_jacobian = state_jacobian
    return stage_finder
end
