using Plots; plotly()
using LinearAlgebra: I, tr
using StatsBase: rmsd
sigdigits = 4

function x_exact_calc(t)
    x = zeros(length(t), 2)
    x[:,1] .= (t.^2 .- t .+ 1.0).*exp.(-t) - exp.(-2.0.*t)
    x[:,2] = (t .+ 1.0).*exp.(-t)
    return x
end

b_calc(t) = [0.0, exp(-t)]
f_calc(t) = [0.5*t*exp(-t), exp(-t)]
dfdt_calc(t) = [0.5*(1.0 - t)*exp(-t), -exp(-t)]
dfdt2_calc(t) = [0.5*(t - 2.0)*exp(-t), exp(-t)]
A_calc(t) = [-2.0 t
             0.0  -1.0]
dAdt_calc(t) = [0.0 1.0
                0.0 0.0]
Ainv_calc(t) = [-0.5 -0.5*t
                0.0  -1.0]
dAinvdt_calc(t) = [0.0 -0.5
                   0.0 0.0]
dAinvdt2_calc(t) = [0.0 0.0
                    0.0 0.0]

function z0_calc(t1, t2)
    # note: first term of Magnus series
    dt = t1 - t2
    # so I need to go to 3rd order to get more accurate (tau-t) expansion
    z0 = [2.0*dt 0.5*(t2^2-t1^2) - dt^3/12.0
          0.0         dt]
    return z0
end

function riemann_sum_right(x0, t_vect)
    nx = length(x0)
    nt = length(t_vect)

    t0 = t_vect[1]
    t_idxs = 2:length(t_vect)
    t_prev = t0

    x = copy(x0)
    x_tmp = copy(x0)
    x_list = [x0...]

    for n in t_idxs
        x_tmp .= x
        t = t_vect[n]
        t_prev = t_vect[n-1]
        dt = t - t_prev

        b1 = b_calc(t)
        z0 = z0_calc(t, t_prev)

        # riemann sum right integration
        x1 = exp(-z0) * x_tmp
        x2 = dt * b1

        x .= x1 + x2
        append!(x_list, x)
    end
    x_list = reshape(x_list, nx, nt) |> transpose
    return x_list
end

function trapezoid_rule(x0, t_vect)
    nx = length(x0)
    nt = length(t_vect)

    t0 = t_vect[1]
    t_idxs = 2:length(t_vect)
    t_prev = t0

    x = copy(x0)
    x_tmp = copy(x0)
    x_list = [x0...]

    for n in t_idxs
        x_tmp .= x
        t = t_vect[n]
        t_prev = t_vect[n-1]
        dt = t - t_prev

        b0 = b_calc(t_prev)
        b1 = b_calc(t)
        z0 = z0_calc(t, t_prev)

        # trapezoid rule integration
        x1 = exp(-z0) * x_tmp
        x2 = 0.5 * dt * (exp(-z0)*b0 + b1)

        x .= x1 + x2
        append!(x_list, x)
    end
    x_list = reshape(x_list, nx, nt) |> transpose
    return x_list
end

function midpoint_rule(x0, t_vect; discrete = true)
    nx = length(x0)
    nt = length(t_vect)

    t0 = t_vect[1]
    t_idxs = 2:length(t_vect)
    t_prev = t0

    x = copy(x0)
    x_tmp = copy(x0)
    x_list = [x0...]

    for n in t_idxs
        if discrete
            t_prev = t_vect[n-1]
            x_tmp .= x
        else
            t_prev = t0
            x_tmp .= x0
        end

        t = t_vect[n]
        dt = t - t_prev
        t_mid = t_prev + dt/2.0

        b = b_calc(t_mid)
        z0 = z0_calc(t, t_prev)
        z_mid = z0_calc(t, t_mid)

        # midpoint rule integration
        # note: more accurate but requires 2 matrix exponentials per step
        x1 = exp(-z0) * x_tmp
        x2 = dt * exp(-z_mid) * b

        x .= x1 + x2
        append!(x_list, x)
    end
    x_list = reshape(x_list, nx, nt) |> transpose
    return x_list
end

function scalar_generator(x0, t_vect; discrete = true, order = 1)
    nx = length(x0)
    nt = length(t_vect)

    t0 = t_vect[1]
    t_idxs = 2:length(t_vect)
    t_prev = t0

    x = copy(x0)
    x_tmp = copy(x0)
    x_list = [x0...]

    for n in t_idxs
        if discrete
            t_prev = t_vect[n-1]
            x_tmp .= x
        else
            t_prev = t0
            x_tmp .= x0
        end
        t = t_vect[n]

        A = A_calc(t)
        Ainv = Ainv_calc(t)
        dAinvdt = dAinvdt_calc(t)
        dAinvdt2 = dAinvdt2_calc(t)
        f = f_calc(t)
        dfdt = dfdt_calc(t)
        dfdt2 = dfdt2_calc(t)
        z0 = z0_calc(t, t_prev)

        xG0 = f + exp(-z0)*(x_tmp - f)

        # dfdt_eff = (I - dAinvdt + 1/nx*tr(dAinvdt*A)*Ainv) \ dfdt   # full recursion
        # dfdt_eff = (I + dAinvdt - 1/nx*tr(dAinvdt*A)*Ainv) * dfdt # truncated correction
        # xG1 = (I - exp(-z0)*(I + 1/nx*tr(z0*Ainv)*A)) * Ainv * dfdt_eff

        # t-expansion is equivalent for this example
        # 1/nx*tr(z0*Ainv) = -(t-t_prev), and tr(dAinvdt*A) = 0
        # dfdt_eff = (I - dAinvdt) \ dfdt
        # xG1 = (I - exp(-z0)*(I - (t-t_prev)*A)) * Ainv * dfdt_eff

        # note: need 2nd Magnus series term to get 3rd order accuracy
        D = (t-t_prev)*exp(-z0) + (I - exp(-z0))*Ainv
        E = -0.5*(t-t_prev)^2*exp(-z0)
        X = inv(I - dAinvdt)
        Y = inv(I - dAinvdt2*X*Ainv - 2.0*dAinvdt)

        C = (E + D*X*Ainv)*Y
        B = (D + C*dAinvdt2)*X

        # 2nd order expansion
        xG1 = B * dfdt
        xG2 = C * dfdt2

        x .= xG0
        order >= 1 ? x .+= xG1 : nothing
        # order >= 2 ? x .+= xG2 : nothing
        x .+= xG2
        append!(x_list, x)
    end
    x_list = reshape(x_list, nx, nt) |> transpose
    return x_list
end

function matrix_generator(x0, t_vect; discrete = true, order = 1)
    nx = length(x0)
    nt = length(t_vect)

    t0 = t_vect[1]
    t_idxs = 2:length(t_vect)
    t_prev = t0

    x = copy(x0)
    x_tmp = copy(x0)
    x_list = [x0...]

    for n in t_idxs
        if discrete
            t_prev = t_vect[n-1]
            x_tmp .= x
        else
            t_prev = t0
            x_tmp .= x0
        end
        t = t_vect[n]

        A = A_calc(t)
        dAdt = dAdt_calc(t)
        Ainv = Ainv_calc(t)
        dAinvdt = dAinvdt_calc(t)
        f = f_calc(t)
        dfdt = dfdt_calc(t)
        dfdt2 = dfdt2_calc(t)
        z0 = z0_calc(t, t_prev)

        xG0 = f + exp(-z0)*(x_tmp - f)
        x .= xG0
        if order == 1
            B = I - exp(-z0)*(I+z0)
            xG1 = B*Ainv*dfdt
            x .+= xG1
        elseif order == 2
            D = I - exp(-z0)*(I+z0)
            E2 = -0.5*exp(-z0)*z0*Ainv*z0
            Q = Ainv*(Ainv*dAdt - dAdt*Ainv)
            X = inv(I - 1.5*Q)
            C = (E2 + D*Ainv)*X
            B = D + 0.5*C*Q*A

            xG1 = B*Ainv*dfdt
            xG2 = C*(dAinvdt*dfdt + Ainv*dfdt2)
            x .+= xG1 .+ xG2
        end
        append!(x_list, x)
    end
    x_list = reshape(x_list, nx, nt) |> transpose
    return x_list
end

dt = 0.5
t0 = 0.0
tf = 5.0
x0 = [0.0, 1.0]
t_vect = t0:dt:tf
x_exact = x_exact_calc(t_vect)
order = 2

x_SG = scalar_generator(x0, t_vect; order)
x_MG = matrix_generator(x0, t_vect; order)
x_RSR = riemann_sum_right(x0, t_vect)
x_TR = trapezoid_rule(x0, t_vect)
x_MR = midpoint_rule(x0, t_vect)
# x = x_SG
x = x_MG
# x = x_RSR
# x = x_TR
# x = x_MR

t_fine = t0:0.01:tf
plt = plot(t_fine, x_exact_calc(t_fine), color = [:black :gray], linewidth = 1.5);
plot!(t_vect, x, color = [:red :blue], linewidth = 1.5, line = :dash)
scatter!(t_vect, x, color = [:red :blue], markersize = [4,4])
display(plt)

println("Root mean squared error:")
println("    scalar generator  = ", round(rmsd(x_SG, x_exact); sigdigits))
println("    matrix generator  = ", round(rmsd(x_MG, x_exact); sigdigits))
println("    riemann sum right = ", round(rmsd(x_RSR, x_exact); sigdigits))
println("    trapezoid rule    = ", round(rmsd(x_TR, x_exact); sigdigits))
println("    midpoint rule     = ", round(rmsd(x_MR, x_exact); sigdigits))

# estimate order of accuracy by halving time step and recalculate RMSEs
x0 = [0.0, 1.0]
t_vect_2 = t0:(dt/2):tf
x_exact_2 = x_exact_calc(t_vect_2)

x_SG_2 = scalar_generator(x0, t_vect_2; order)
x_MG_2 = matrix_generator(x0, t_vect_2; order)
x_RSR_2 = riemann_sum_right(x0, t_vect_2)
x_TR_2 = trapezoid_rule(x0, t_vect_2)
x_MR_2 = midpoint_rule(x0, t_vect_2)

p_SG = log2(rmsd(x_SG, x_exact) / rmsd(x_SG_2, x_exact_2))
p_MG = log2(rmsd(x_MG, x_exact) / rmsd(x_MG_2, x_exact_2))
p_RSR = log2(rmsd(x_RSR, x_exact) / rmsd(x_RSR_2, x_exact_2))
p_TR = log2(rmsd(x_TR, x_exact) / rmsd(x_TR_2, x_exact_2))
p_MR = log2(rmsd(x_MR, x_exact) / rmsd(x_MR_2, x_exact_2))

println("\nOrder of accuracy:")
println("    scalar generator  = ", round(p_SG; sigdigits))
println("    matrix generator  = ", round(p_MG; sigdigits))
println("    riemann sum right = ", round(p_RSR; sigdigits))
println("    trapezoid rule    = ", round(p_TR; sigdigits))
println("    midpoint rule     = ", round(p_MR; sigdigits))
println("\ndone")
