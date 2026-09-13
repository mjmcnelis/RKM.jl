using OrdinaryDiffEq, LambertW
using OrdinaryDiffEqBDF
using StatsBase: rmsd, mean
using ForwardDiff: Dual, partials
using Plots; plotly(); #gr()

show_plot = false

# logistic
# f(y, t) = (y + 0.5)*(0.5 - y)
# J_y(y, t) = -2.0 * y
# function f_ode(du, u, p, t)
#     du[1] = (u[1] + 0.5)*(0.5 - u[1])
#     return nothing
# end
# y_exact(y0, t_vect) = @. 0.5 * tanh(t_vect/2 + atanh(2*y0))

# flame ODE
f(y, t) = y^2 * (1.0 - y)
J_y(y, t) = y * (2.0 - 3.0*y)
function f_ode(du, u, p, t)
    du[1] = u[1]^2 * (1.0 - u[1])
    return nothing
end
y_exact(y0, t_vect) = @. 1.0 / (1.0 + lambertw(((1.0/y0) - 1.0) * exp(((1.0/y0) - 1.0) - t_vect)))

function magnus_step(y_n, t_n, dt, y_prev, n; tol = 1e-8, maxiter = 50)
    t = t_n + dt

    f_n = f(y_n, t_n)
    J_n = J_y(y_n, t_n)

    # initial guess (tangent quadratic)
    y = y_n + f_n*dt + J_n*f_n*dt^2/2.0

    trajectory = "quadratic"
    # trajectory = "linear"

    # tangent = true
    tangent = false

    iter = 0
    for k = 1:maxiter
        # reconstruct trajectory:
        #   y(s) = y_n + b*s + a*s^2
        if trajectory == "quadratic"
            if n == 1 || tangent
                # quadratic curve tangent to f_n at t_n
                b = f_n
                a = (y - y_n - dt*f_n)/dt^2
            else
                # quadratic interpolation through (y_prev, y_n, y)
                b = (y - y_prev)/(2.0*dt)
                a = (y - 2.0*y_n + y_prev) / (2.0*dt^2)
            end
        elseif trajectory == "linear"
            # linear interpolation through (y_n, y)
            b = (y - y_n) / dt
            a = 0.0
        end

        # trajectory derivatives at end point
        dydt = b + 2.0*a*dt
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

        # eta = D0*g
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
        n_prev = max(1, n-1)
        y[n+1], iters[n] = magnus_step(y[n], t[n], dt, y[n_prev], n)
    end
    return t, y, iters
end

# y0 = 0.0
# t0 = 0.0
# tf = 10.0
# dt = 0.01

δ = 1e-4
y0 = δ
t0 = 0.0
tf = 2.0/δ
dt = 0.01

@time t, y_mag, iters = solve_magnus(y0, t0, tf, dt)
GC.gc()

prob = ODEProblem(f_ode, [y0], (t0, tf))
# method = Trapezoid
method = QBDF2
@time sol = solve(prob, method(); maxiters = 10^7, adaptive = false, dt, saveat = t)
y_tr = hcat(sol.u...)'

y_ex = y_exact(y0, t)

println("\nmax iterations = ", maximum(iters))
println("average iterations = ", mean(iters), "\n")
println("Root mean squared error:")
println("  $method = ", round(rmsd(y_tr, y_ex), sigdigits = 4),)
println("  Magnus = ", round(rmsd(y_mag, y_ex), sigdigits = 4),)

if show_plot
    plt = plot(t, y_ex, label = "Exact", lw = 2, legend = :outertopright, size = (900,500),);
    plot!(sol.t, y_tr, #=marker = :circle,=# label = "OrdinaryDiffEq Trapezoid");
    plot!(t, y_mag, #=marker = :circle,=# linestyle = :dash, label = "Magnus iteration");
    display(plt)
end
GC.gc()
println("\ndone")