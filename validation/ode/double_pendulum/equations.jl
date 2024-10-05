
const g = 9.81

# for RKM
function dy_dt!(f, y, t; p, kwargs...)
    L = p[1]
    θ1, θ2, pθ1, pθ2 = y

    cos12 = cos(θ1 - θ2)
    factor = 6.0/(16.0 - 9.0*cos12^2)

    # d/dt(theta_1)
    f[1] = (2.0*pθ1 - 3.0*pθ2*cos12) * factor
    # d/dt(theta_2)
    f[2] = (8.0*pθ2 - 3.0*pθ1*cos12) * factor

    θdot1 = f[1]
    θdot2 = f[2]
    θdot1θdot2sin12 = θdot1 * θdot2 * sin(θ1 - θ2)

    # note: normalized canonical momentum by moment of inertia
    # d/dt(p_theta_1) / (m*L^2)
    f[3] = -(θdot1θdot2sin12 + 3.0*sin(θ1)*g/L) / 2.0
    # d/dt(p_theta_2) / (m*L^2)
    f[4] = (θdot1θdot2sin12 - sin(θ2)*g/L) / 2.0
    return nothing
end

# for OrdinaryDiffEq
function dy_dt!(f, y, p, t)
    L = p[1]
    θ1, θ2, pθ1, pθ2 = y

    cos12 = cos(θ1 - θ2)
    factor = 6.0/(16.0 - 9.0*cos12^2)

    # d/dt(theta_1)
    f[1] = (2.0*pθ1 - 3.0*pθ2*cos12) * factor
    # d/dt(theta_2)
    f[2] = (8.0*pθ2 - 3.0*pθ1*cos12) * factor

    θdot1 = f[1]
    θdot2 = f[2]
    θdot1θdot2sin12 = θdot1 * θdot2 * sin(θ1 - θ2)

    # note: normalized canonical momentum by moment of inertia
    # d/dt(p_theta_1) / (m*L^2)
    f[3] = -(θdot1θdot2sin12 + 3.0*sin(θ1)*g/L) / 2.0
    # d/dt(p_theta_2) / (m*L^2)
    f[4] = (θdot1θdot2sin12 - sin(θ2)*g/L) / 2.0
    return nothing
end
