
@kwdef struct PIDControlBeta{T <: AbstractFloat}
    # why not just always use Float64 for betas?
    beta1::T
    beta2::T
    beta3::T
    alpha2::T = 0.0
    alpha3::T = 0.0
end

function PIDControlK(; kI, kP, kD, alpha2 = 0.0, alpha3 = 0.0)
    # relation between beta and k control parameters
    beta1 = kI + kP + kD
    beta2 = -kP - 2kD
    beta3 = kD
    return PIDControlBeta(; beta1, beta2, beta3, alpha2, alpha3)
end

function BasicControl()
    return PIDControlK(; kI = 1.0, kP = 0.0, kD = 0.0, alpha2 = 0.0, alpha3 = 0.0)
end

function PIControl(; kI = 0.3, kP = 0.4)
    return PIDControlK(; kI, kP, kD = 0.0, alpha2 = 0.0, alpha3 = 0.0)
end

function H312Control(; kI = 2/9)
    beta1 = kI/4
    beta2 = kI/2
    beta3 = kI/4
    return PIDControlBeta(; beta1, beta2, beta3, alpha2 = 0.0, alpha3 = 0.0)
end

function H321PredictiveControl(; kI = 0.1, kP = 0.45)
    beta1 = 3kI/4 + kP/2
    beta2 = kI/2
    beta3 = -kI/4 - kP/2
    return PIDControlBeta(; beta1, beta2, beta3, alpha2 = -1.0, alpha3 = 0.0)
end

function H211Control(; kI = 1/6)
    beta1 = kI/2
    beta2 = kI/2
    return PIDControlBeta(; beta1, beta2, beta3 = 0.0, alpha2 = 0.0, alpha3 = 0.0)
end

function H211bPredictiveControl(; b = 4.0)
    2 <= b <= 8 ? nothing : @warn "Recommend using b ∈ [2,8]"
    beta1 = 1/b
    beta2 = 1/b
    alpha2 = 1/b
    # TODO: use alpha1 and alpha2 instead of predictive = true
    return PIDControlBeta(; beta1, beta2, beta3 = 0.0, alpha2, alpha3 = 0.0)
end