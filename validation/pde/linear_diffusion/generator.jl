using Revise, RKM, LinearSolve, ForwardDiff
using AppleAccelerate
using Plots; plotly()
if !(@isdefined dy_dt!)
    function dy_dt!(f, y, t; p, abstract_params)
        # note: here dx is not a sensitivity parameter
        a = p[1]
        dx = abstract_params
        N = length(y)
        for i in 1:N
            # Dirichlet boundary conditions: y[0] = 0.0, y[N+1] = 0.0
            ym = i == 1 ? 0.0 : y[i-1]
            yp = i == N ? 0.0 : y[i+1]
            yc = y[i]
            f[i] = (yp - 2.0*yc + ym)*a/dx^2
        end
        nothing
    end
    function gauss(x, t; p, t0)
        a = p[1]
        return sqrt(t0/t) * exp(-x^2/(4.0*a*t))
    end
end

# show_plot = true
show_plot = false

a = 0.25                    # diffusion constant
Nx = 21
x = range(-10, 10, Nx)      # grid points
dx = x[2] - x[1]            # uniform spacing
p = [a]                     # parameters
abstract_params = dx        # non-sensitivity parameters

t0 = 1.0                    # initial conditions
y0 = gauss.(x, t0; p, t0)

tf = 100.0
dt0 = 0.01

CFL = 2.0*a*dt0/dx^2        # CFL number
@show CFL
Nt = 300                    # temporal stride (for plot)

# generate sparsity pattern via nansafe
sparsity = nansafe_state_jacobian(y0, t0, dy_dt!, p; chunk_size = 1, abstract_params)
display(sparsity)

# test_nansafe(; x = ForwardDiff.Dual(NaN, 1.0))
# test_nansafe(; x = ForwardDiff.Dual(NaN, 0.0))
# q()

options = SolverOptions(
              method = BackwardEuler1(),
              adaptive = Fixed(),
              state_jacobian = FiniteJacobian(),
              sensitivity = DecoupledDirect(),
          )

@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p; abstract_params)

_, y = get_solution(sol)
_, S = get_sensitivity(sol)

using LinearAlgebra: I, inv, det, lu, lu!, ldiv!, eigvals

A = Matrix(sparsity) #- 1e-3I   # even if regulator can reduce noise, the curves shifts wrong
G(t, t0) = exp(A*(t - t0))
z(t, t0) = -A*(t - t0)

(; nt, ny, np) = get_dimensions(sol)

y = y0 |> deepcopy
ymm = y0 |> deepcopy
ym = y0 |> deepcopy
yp = y0 |> deepcopy
ypp = y0 |> deepcopy
f = zeros(ny)
dfdp = zeros(ny, 1)
dfdpm = zeros(ny, 1)
dfdpmm = zeros(ny, 1)
dfdpp = zeros(ny, 1)
dfdppp = zeros(ny, 1)
ddtdfdp = zeros(ny, 1)
d2dt2dfdp = zeros(ny, 1)

ode_wrap_p! = RKM.ODEWrapperParam([t0], y, abstract_params, dy_dt!)

sensitivity = options.sensitivity
sensitivity = RKM.reconstruct_sensitivity(sensitivity, ode_wrap_p!, f, p, false)
param_jacobian = sensitivity.param_jacobian

dfdp_list = Float64[]
dydp_list = Float64[]

dn = 1

skip = 10                   # skipping time points is a big advantage
# t_idxs = 2:skip:nt-1
# t_idxs = dn+1:skip:nt-dn
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
    # f not used anymore... (revisit old jacobian methods)
    RKM.evaluate_jacobian!(param_jacobian, dfdp, ode_wrap_p!, p)
    # minus-minus
    RKM.set_wrapper!(ode_wrap_p!, tmm, ymm)
    RKM.evaluate_jacobian!(param_jacobian, dfdpmm, ode_wrap_p!, p)
    # minus
    RKM.set_wrapper!(ode_wrap_p!, tm, ym)
    RKM.evaluate_jacobian!(param_jacobian, dfdpm, ode_wrap_p!, p)
    # plus
    RKM.set_wrapper!(ode_wrap_p!, tp, yp)
    RKM.evaluate_jacobian!(param_jacobian, dfdpp, ode_wrap_p!, p)
    # plus-plus
    RKM.set_wrapper!(ode_wrap_p!, tpp, ypp)
    RKM.evaluate_jacobian!(param_jacobian, dfdppp, ode_wrap_p!, p)

    # 1st time derivative
    # ddtdfdp .= (dfdp - dfdpm) / dt             # backward
    # ddtdfdp .= (dfdpp - dfdp) / dt             # forward
    @. ddtdfdp = dfdpp - dfdpm                   # central
    @. ddtdfdp /= (2.0*dt)
    # @. ddtdfdp = -dfdppp + 8.0*dfdpp - 8.0*dfdpm + dfdpmm
    # @. ddtdfdp /= (12.0*dt)

    # 2nd time derivative (classic)
    @. d2dt2dfdp = dfdpp - 2.0*dfdp + dfdpm
    @. d2dt2dfdp /= dt^2

    # 2nd time derivative (5-point), this was more noisy...
    # @. d2dt2dfdp = -dfdppp + 16.0*dfdpp - 30.0*dfdp + 16.0*dfdpm - dfdpmm
    # @. d2dt2dfdp /= (12.0*dt^2)

    append!(dfdp_list, dfdp)

    # TODO: fastest way to determine if sparse, singular A, get min eigenvalue?

    G0 = G(t, t0)
    z0 = z(t, t0)

    #=
    Gamma1 = I - G0
    Gamma2 = I - G0*(I + z0)
    Gamma3 = I - G0*(I + z0 + 0.5*z0*z0)

    dydp0 = -Gamma1 * (A \ dfdp)
    dydp1 = -Gamma2 * (A \ (A \ ddtdfdp))
    dydp2 = -Gamma3 * (A \ (A \ (A \ d2dt2dfdp)))
    =#

    F = lu(A)             # note: don't use lu! unless reset A

    S0 = -dfdp
    ldiv!(F, S0)

    dS1 = -ddtdfdp
    ldiv!(F, dS1)
    ldiv!(F, dS1)

    dS2 = -d2dt2dfdp
    ldiv!(F, dS2)
    ldiv!(F, dS2)
    ldiv!(F, dS2)

    # current favorite b/c reduces projections
    dydp0 = S0 + dS1 + dS2
    dydp1 = -G0*(S0 + dS1 + dS2 + z0*(dS1 + dS2 + 0.5*z0*dS2))
    # dydp2 = 0*dfdp

    # converging, but 2nd-order term pretty noisy
    dydp = dydp0 + dydp1 #+ dydp2

    append!(dydp_list, dydp)
end

dfdp_list = reshape(dfdp_list, ny, length(t_idxs)) |> transpose
dydp_list = reshape(dydp_list, ny*np, length(t_idxs)) |> transpose

# plot(sol.t[2:end-1], dfdp_list) |> display

plt = plot(sol.t, S, legend = :outertopright, size = (900,500),);
plot!(sol.t[t_idxs], dydp_list, label = "", color = :black, line = :dash);
display(plt)

GC.gc()
println("\ndone")
