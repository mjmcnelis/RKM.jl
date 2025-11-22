using Revise, RKM, LinearSolve, ForwardDiff
using AppleAccelerate
using DoubleFloats
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/pde/linear_diffusion/equations_dirichlet.jl") : nothing

precision = Float64
# precision = Double64

a = 0.25                    # diffusion constant
Nx = 21
x = range(-10, 10, Nx)      # grid points
dx = x[2] - x[1]            # uniform spacing
p = [a]                     # parameters
abstract_params = dx        # non-sensitivity parameters
# p = [a, dx]               # works for both sensitivity parameters
# abstract_params = nothing

t0 = 1.0                    # initial conditions
y0 = gauss.(x, t0; p, t0)

tf = 100.0
dt0 = 0.01

CFL = 2.0*a*dt0/dx^2        # CFL number
@show CFL
Nt = 300                    # temporal stride (for plot)

# generate sparsity pattern via nansafe
sparsity = nansafe_state_jacobian(y0, t0, dy_dt!, p; chunk_size = 1, abstract_params)

if precision == Float64
    method = BackwardEuler1()
else
    method = Ketcheson4()
end

options = SolverOptions(;
              method, adaptive = Fixed(),
              state_jacobian = FiniteJacobian(),
              sensitivity = DecoupledDirect(; param_jacobian = ForwardJacobian(),),
              precision,)

@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p; abstract_params)
_, y = get_solution(sol)
_, S = get_sensitivity(sol)

using LinearAlgebra: I, mul!, lu, lu!, ldiv!, eigvals

A = Matrix(sparsity)
G(t, t0) = exp(A*(t - t0))
z(t, t0) = -A*(t - t0)

(; nt, ny, np) = get_dimensions(sol)

y = y0 .|> precision
ymm = y |> deepcopy
ym = y |> deepcopy
yp = y |> deepcopy
ypp = y |> deepcopy
dydp = zeros(precision, ny, np)
dSG = zeros(precision, ny, np)

f = zeros(precision, ny)
dfdp = zeros(precision, ny, np)
dfdpm = zeros(precision, ny, np)
dfdpmm = zeros(precision, ny, np)
dfdpp = zeros(precision, ny, np)
dfdppp = zeros(precision, ny, np)

ddtdfdp = zeros(precision, ny, np)
d2dt2dfdp = zeros(precision, ny, np)
d3dt3dfdp = zeros(precision, ny, np)
d4dt4dfdp = zeros(precision, ny, np)

t0 = precision(t0)
p = precision.(p)

ode_wrap_p! = RKM.ODEWrapperParam([t0], y, abstract_params, dy_dt!)

sensitivity = options.sensitivity
sensitivity = RKM.reconstruct_sensitivity(sensitivity, ode_wrap_p!, f, p, false)
param_jacobian = sensitivity.param_jacobian

dydp_list = Float64[]

dn = 10                      # stride used for finite difference
skip = 100                   # skipping time points is a big advantage
t_idxs = 2dn+1:skip:nt-2dn

@time for n in t_idxs
    tmm = sol.t[n-2dn]
    tm = sol.t[n-dn]
    t = sol.t[n]
    tp = sol.t[n+dn]
    tpp = sol.t[n+2dn]

    dt = tp - t

    ymm .= view(sol.y, (1 + (n-1-2dn)*ny):(n-2dn)*ny)
    ym .= view(sol.y, (1 + (n-1-dn)*ny):(n-dn)*ny)
    y .= view(sol.y, (1 + (n-1)*ny):n*ny)
    yp .= view(sol.y, (1 + (n-1+dn)*ny):(n+dn)*ny)
    ypp .= view(sol.y, (1 + (n-1+2dn)*ny):(n+2dn)*ny)

    # current
    RKM.set_wrapper!(ode_wrap_p!, t, y)
    RKM.evaluate_jacobian!(param_jacobian, dfdp, ode_wrap_p!, p, f)
    # minus-minus
    RKM.set_wrapper!(ode_wrap_p!, tmm, ymm)
    RKM.evaluate_jacobian!(param_jacobian, dfdpmm, ode_wrap_p!, p, f)
    # minus
    RKM.set_wrapper!(ode_wrap_p!, tm, ym)
    RKM.evaluate_jacobian!(param_jacobian, dfdpm, ode_wrap_p!, p, f)
    # plus
    RKM.set_wrapper!(ode_wrap_p!, tp, yp)
    RKM.evaluate_jacobian!(param_jacobian, dfdpp, ode_wrap_p!, p, f)
    # plus-plus
    RKM.set_wrapper!(ode_wrap_p!, tpp, ypp)
    RKM.evaluate_jacobian!(param_jacobian, dfdppp, ode_wrap_p!, p, f)

    # 0.012555 seconds (25.17 k allocations: 5.926 MiB, 54.70% compilation time)

    # time derivatives (central differences)
    @. ddtdfdp = dfdpp - dfdpm                                          # 1st derivative
    @. ddtdfdp /= (2.0*dt)

    @. d2dt2dfdp = dfdpp - 2.0*dfdp + dfdpm                             # 2nd derivative
    @. d2dt2dfdp /= dt^2

    @. d3dt3dfdp = dfdppp - 2.0*dfdpp + 2.0*dfdpm - dfdpmm              # 3rdd derivative
    @. d3dt3dfdp /= (2.0*dt^3)

    @. d4dt4dfdp = dfdppp - 4.0*dfdpp + 6.0*dfdp - 4.0*dfdpm + dfdpmm   # 4th derivative
    @. d4dt4dfdp /= dt^4

    G0 = G(t, t0)
    z0 = z(t, t0)

    F = lu(A)             # note: don't use lu! unless reset A

    S0 = dfdp
    S0 .*= -1.0
    ldiv!(F, S0)

    dS1 = ddtdfdp
    dS1 .*= -1.0
    ldiv!(F, dS1)
    ldiv!(F, dS1)

    dS2 = d2dt2dfdp
    dS2 .*= -1.0
    ldiv!(F, dS2)
    ldiv!(F, dS2)
    ldiv!(F, dS2)

    dS3 = d3dt3dfdp
    dS3 .*= -1.0
    ldiv!(F, dS3)
    ldiv!(F, dS3)
    ldiv!(F, dS3)
    ldiv!(F, dS3)

    dS4 = d4dt4dfdp
    dS4 .*= -1.0
    ldiv!(F, dS4)
    ldiv!(F, dS4)
    ldiv!(F, dS4)
    ldiv!(F, dS4)
    ldiv!(F, dS4)

    # TODO: fastest way to determine if sparse, singular A, get min eigenvalue?

    # is there a corresponding ODE for Gamma that I can project onto?
    # maybe its form is simpler to deal with
    #=
    Gamma1 = I - G0
    Gamma2 = I - G0*(I + z0)
    Gamma3 = I - G0*(I + z0 + z0^2/2.0)
    Gamma4 = I - G0*(I + z0 + z0^2/2.0 + z0^3/6.0)
    Gamma5 = I - G0*(I + z0 + z0^2/2.0 + z0^3/6.0 + z0^4/24.0)

    dydp .= 0.0
    mul!(dSG, Gamma1, S0)
    dydp .+= dSG
    mul!(dSG, Gamma2, dS1)
    dydp .+= dSG
    mul!(dSG, Gamma3, dS2)
    dydp .+= dSG
    mul!(dSG, Gamma4, dS3)
    dydp .+= dSG
    mul!(dSG, Gamma5, dS4)
    dydp .+= dSG
    =#

    # not sure if I should keep going with this
    # z0dS1 = (t - t0) * (A \ ddtdfdp)
    # z0dS2 = (t - t0) * (A \ (A \ d2dt2dfdp))
    # z02dS2 = -(t - t0)^2 * (A \ d2dt2dfdp)
    # dydp0 = (I - G0) * (S0 + dS1 + dS2)
    # dydp1 = -G0*(z0dS1 + z0dS2 + 0.5*z02dS2)
    # dydp2 = 0*dfdp

    # dydp0 = -(I - G0) * (A \ (dfdp + A \ (ddtdfdp + A \ d2dt2dfdp)))
    # dydp1 = -(t - t0) * G0 * (A \ (ddtdfdp + A \ d2dt2dfdp))
    # dydp2 = 0.5 * (t - t0)^2 * G0 * (A \ d2dt2dfdp)

    # can you reuse a matrix factorization if A is sparse?

    # current favorite b/c reduces projections
    # dydp0 = S0 + dS1 + dS2
    # dydp1 = -G0*(S0 + dS1 + dS2 + z0*(dS1 + dS2 + 0.5*z0*dS2))
    # dydp2 = 0*dfdp

    # besides reducing projections, can you further
    # group dS terms to reduce linear solves?

    dydp .= 0.0
    @. dydp += S0 + dS1 + dS2 + dS3 + dS4
    dydp1 = -G0*(S0 + dS1 + dS2 + dS3 + dS4 +
                 z0*(dS1 + dS2 + dS3 + dS4 +
                     z0*((dS2 + dS3 + dS4)/2.0 +
                          z0*((dS3 + dS4)/6.0 +
                               z0*dS4/24.0
                             )
                        )
                    )
                )
    dydp .+= dydp1

    append!(dydp_list, dydp)
end

dydp_list = reshape(dydp_list, ny*np, length(t_idxs)) |> transpose

plt = plot(sol.t, S, legend = :outertopright, size = (900,500),);
plot!(sol.t[t_idxs], dydp_list, label = "", color = :black, line = :dash);
display(plt)

GC.gc()
println("\ndone")
