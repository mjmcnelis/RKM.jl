
abstract type StageFinder end

# good enough start (wrap caches later)
@kwdef struct ImplicitStageFinder{JM, AF, EMM} <: StageFinder where {
                                                               JM <: JacobianMethod,
                                                               AF <: AbstractFactorization,
                                                               EMM <: EigenMaxMethod}
    state_jacobian::JM = FiniteJacobian()
    linear_method::AF = LUFactorization()
    eigenmax::EMM = NoEigenMax()
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
