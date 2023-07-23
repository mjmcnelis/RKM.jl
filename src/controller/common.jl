
function BasicControl(; safety = 0.8)
    return PIDControlK(; kI = 1.0, kP = 0.0, kD = 0.0, predictive = false)
end

function PI34Control(; safety = 0.8)
    return PIDControlK(; kI = 0.3, kP = 0.4, kD = 0.0, predictive = false)
end

function PI33Control(; safety = 0.8)
    return PIDControlK(; kI = 0.3, kP = 0.3, kD = 0.0, predictive = false)
end

function PI42Control(; safety = 0.8)
    return PIDControlK(; kI = 0.4, kP = 0.2, kD = 0.0, predictive = false)
end

function H312Control(; kI = 2/9, safety = 0.8)
    beta1 = kI/4
    beta2 = kI/2
    beta3 = kI/4
    return PIDControlBeta(; beta1, beta2, beta3, predictive = false)
end

function H321PredictiveControl(; kI = 0.1, kP = 0.45, safety = 0.8)
    beta1 = 3kI/4 + kP/2
    beta2 = kI/2
    beta3 = -kI/4 - kP/2
    return PIDControlBeta(; beta1, beta2, beta3, predictive = true)
end