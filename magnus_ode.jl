using OrdinaryDiffEq, LambertW
using StatsBase: rmsd, mean
using ForwardDiff: Dual, partials
using Plots; plotly()
# gr()

# logistic
f(y, t) = (y + 0.5)*(0.5 - y)
J_y(y, t) = -2.0 * y
# OrdinaryDiffEq
function f_ode(du, u, p, t)
    du[1] = (u[1] + 0.5)*(0.5 - u[1])
    return nothing
end

# flame ODE
# f(y, t) = y^2 - y^3
# J_y(y, t) = 2.0*y - 3.0*y^2

function magnus_step(y_n, t_n, dt; tol = 1e-8, maxiter = 50)
    f_n = f(y_n, t_n)
    J_n = J_y(y_n, t_n)

    # initial guess (tangent quadratic)
    y = y_n + f_n*dt + J_n*f_n*dt^2/2.0
    t = t_n + dt

    iter = 0
    for k = 1:maxiter
        # reconstruct quadratic trajectory
        # y(s) = y_n + f_n*s + a*s^2

        # trajectory derivatives at end point
        a = (y - y_n - dt*f_n)/dt^2
        dydt = f_n + 2.0*a*dt
        dydt2 = 2.0*a

        # Jacobian and residual source term
        J = J_y(y, t)
        r = f(y, t) - dydt

        y_dual = Dual(y, dydt)
        t_dual = Dual(t, 1.0)

        f_dual = f(y_dual, t_dual)
        J_dual = J_y(y_dual, t_dual)

        dJdt = partials(J_dual, 1)
        drdt = partials(f_dual, 1) - dydt2

        g = -J \ r
        dgdt = -J \ (dJdt*g + drdt)

        # scalar propagator z0 = -∫J(t)dt (approx w/ trapezoid rule)
        z0 = -(J_n + J)*dt/2.0

        # matrix generator formula
        D0 = -expm1(-z0)
        D1 = D0 - exp(-z0)*z0
        eta = D0*g + D1*(J \ dgdt)

        y += eta
        iter += 1
        if abs(eta) < tol*(1.0 + abs(y))
            break
        end
    end

    return y, iter
end

function solve_magnus(y0, t0, tf, dt)
    t = collect(t0:dt:tf)

    y = zeros(length(t))
    iters = zeros(Int, length(t)-1)

    y[1] = y0
    for n = 1:length(t)-1
        y[n+1], iters[n] = magnus_step(y[n], t[n], dt)
    end
    return t, y, iters
end

# ------------------------------------------------------------
# Test
# ------------------------------------------------------------

y0 = 0.0
t0 = 0.0
tf = 10.0
dt = 0.01

# δ = 1e-4
# y0 = δ
# t0 = 0.0
# tf = 2.0/δ
# # dt = δ
# dt = 1.0

@time t, y_mag, iters = solve_magnus(y0, t0, tf, dt)
GC.gc()

prob = ODEProblem(f_ode, [y0], (t0, tf))
@time sol = solve(prob, Trapezoid(); maxiters = 10^7, adaptive = false, dt, saveat = t)
y_tr = hcat(sol.u...)'

# Exact solution
y_exact = @. 0.5*tanh(t/2.0)
# y_exact = @. 1.0 / (1.0 + lambertw(((1.0/y0) - 1.0) * exp(((1.0/y0) - 1.0) - t)))

plt = plot(t, y_exact, label = "Exact", lw = 2, legend = :outertopright, size = (900,500),);
plot!(sol.t, y_tr, #=marker = :circle,=# label = "OrdinaryDiffEq Trapezoid");
plot!(t, y_mag, #=marker = :circle,=# linestyle = :dash, label = "Magnus iteration");
display(plt)

# plt = plot(t, y_exact-y_tr, label = "exact - trapezoid", lw = 1.5, legend = :outertopright, size = (900,500),);
# plot!(t, y_exact-y_mag, label = "exact - magnus", lw = 1.5, linestyle = :dash);
# display(plt)

println("\naverage iterations = ", mean(iters), "\n")

println("Root mean squared error:")
println("    Trapezoid = ", round(rmsd(y_tr, y_exact), sigdigits = 4),)
println("    Magnus    = ", round(rmsd(y_mag, y_exact), sigdigits = 4),)

GC.gc()
println("\ndone")