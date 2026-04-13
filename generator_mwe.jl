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
Ainv_calc(t) = [-0.5 -0.5*t
                0.0  -1.0]
dAinvdt_calc(t) = [0.0 -0.5
                   0.0 0.0]

function z0_calc(t1, t2)
    # note: first term of Magnus series
    z0 = [2.0*(t1-t2) 0.5*(t2^2-t1^2)
          0.0         t1-t2]
    return z0
end

function riemann_sum_left(x0, t_vect)
    nx = length(x0)
    nt = length(t_vect)

    t0 = t_vect[1]
    t_idxs = 2:length(t_vect)
    t_prev = t0

    x = copy(x0)
    x_tmp = copy(x0)
    xRL = [x0...]

    for n in t_idxs
        x_tmp .= x
        t = t_vect[n]
        t_prev = t_vect[n-1]
        dt = t - t_prev

        b0 = b_calc(t_prev)
        z0 = z0_calc(t, t_prev)

        # riemann sum left integration
        x1 = exp(-z0) * x_tmp
        x2 = dt * exp(-z0) * b0

        x .= x1 + x2
        append!(xRL, x)
    end
    xRL = reshape(xRL, nx, nt) |> transpose
    return xRL
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

function affine_midpoint(x0, t_vect)
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
        t_prev = t_vect[n-1]
        t = t_vect[n]
        dt = t - t_prev
        t_mid = t_prev + dt/2.0

        A = A_calc(t_mid)
        Ainv = Ainv_calc(t_mid)
        b = b_calc(t_mid)

        x1 = exp(dt*A) * x_tmp
        x2 = (exp(dt*A) - I) * Ainv * b

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
        f = f_calc(t)
        dfdt = dfdt_calc(t)
        z0 = z0_calc(t, t_prev)

        xG0 = f + exp(-z0)*(x_tmp - f)

        dfdt_eff = (I - dAinvdt + 1/nx*tr(dAinvdt*A)*Ainv) \ dfdt
        xG1 = (I - exp(-z0)*(I + 1/nx*tr(z0*Ainv)*A)) * Ainv * dfdt_eff

        x .= xG0
        order >= 1 ? x .+= xG1 : nothing
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
        Ainv = Ainv_calc(t)
        dAinvdt = dAinvdt_calc(t)
        f = f_calc(t)
        dfdt = dfdt_calc(t)
        dfdt2 = dfdt2_calc(t)
        z0 = z0_calc(t, t_prev)

        xG0 = f + exp(-z0)*(x_tmp - f)
        xG1 = (I - exp(-z0)*(I + z0))*Ainv*dfdt

        dAinvdt_eff = dAinvdt + 0.5*Ainv*(A*dAinvdt - dAinvdt*A)
        xG2 = (I - exp(-z0)*(I + z0 + 0.5*z0^2))*Ainv*(dAinvdt_eff*dfdt + Ainv*dfdt2)

        # this was worse
        # dAinvdt_eff = dAinvdt + 0.5*Ainv*(A*dAinvdt - dAinvdt*A)
        # C = dAinvdt_eff*dfdt + Ainv*dfdt2
        # xG2 = Ainv*C - exp(-z0)*(Ainv +
        #                          0.5*z0*Ainv + 0.5*Ainv*z0 +
        #                          0.5*z0*Ainv*z0)*C

        x .= xG0
        order >= 1 ? x .+= xG1 : nothing
        order >= 2 ? x .+= xG2 : nothing
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

x_SG = scalar_generator(x0, t_vect; order = 1)
x_MG = matrix_generator(x0, t_vect; order = 1)
x_RSL = riemann_sum_left(x0, t_vect)
x_RSR = riemann_sum_right(x0, t_vect)
x_TR = trapezoid_rule(x0, t_vect)
x_MR = midpoint_rule(x0, t_vect)
x_AM = affine_midpoint(x0, t_vect)
# x = x_SG
x = x_MG
# x = x_RSL
# x = x_RSR
# x = x_TR
# x = x_MR
# x = x_AM

t_fine = t0:0.01:tf
plt = plot(t_fine, x_exact_calc(t_fine), color = [:black :gray], linewidth = 1.5);
plot!(t_vect, x, color = [:red :blue], linewidth = 1.5, line = :dash)
scatter!(t_vect, x, color = [:red :blue], markersize = [4,4])
display(plt)

println("Root mean squared error:")
println("    scalar generator  = ", round(rmsd(x_SG, x_exact); sigdigits))
println("    matrix generator  = ", round(rmsd(x_MG, x_exact); sigdigits))
println("    riemann sum left  = ", round(rmsd(x_RSL, x_exact); sigdigits))
println("    riemann sum right = ", round(rmsd(x_RSR, x_exact); sigdigits))
println("    trapezoid rule    = ", round(rmsd(x_TR, x_exact); sigdigits))
println("    midpoint rule     = ", round(rmsd(x_MR, x_exact); sigdigits))
println("    affine midpoint   = ", round(rmsd(x_AM, x_exact); sigdigits))

# estimate order of accuracy by halving time step and recalculate RMSEs
x0 = [0.0, 1.0]
t_vect_2 = t0:(dt/2):tf
x_exact_2 = x_exact_calc(t_vect_2)

x_SG_2 = scalar_generator(x0, t_vect_2; order = 1)
x_MG_2 = matrix_generator(x0, t_vect_2; order = 1)
x_RSL_2 = riemann_sum_left(x0, t_vect_2)
x_RSR_2 = riemann_sum_right(x0, t_vect_2)
x_TR_2 = trapezoid_rule(x0, t_vect_2)
x_MR_2 = midpoint_rule(x0, t_vect_2)
x_AM_2 = affine_midpoint(x0, t_vect_2)

p_SG = log2(rmsd(x_SG, x_exact) / rmsd(x_SG_2, x_exact_2))
p_MG = log2(rmsd(x_MG, x_exact) / rmsd(x_MG_2, x_exact_2))
p_RSL = log2(rmsd(x_RSL, x_exact) / rmsd(x_RSL_2, x_exact_2))
p_RSR = log2(rmsd(x_RSR, x_exact) / rmsd(x_RSR_2, x_exact_2))
p_TR = log2(rmsd(x_TR, x_exact) / rmsd(x_TR_2, x_exact_2))
p_MR = log2(rmsd(x_MR, x_exact) / rmsd(x_MR_2, x_exact_2))
p_AM = log2(rmsd(x_AM, x_exact) / rmsd(x_AM_2, x_exact_2))

println("\nOrder of accuracy:")
println("    scalar generator  = ", round(p_SG; sigdigits))
println("    matrix generator  = ", round(p_MG; sigdigits))
println("    riemann sum left  = ", round(p_RSL; sigdigits))
println("    riemann sum right = ", round(p_RSR; sigdigits))
println("    trapezoid rule    = ", round(p_TR; sigdigits))
println("    midpoint rule     = ", round(p_MR; sigdigits))
println("    affine midpoint   = ", round(p_AM; sigdigits))
println("\ndone")
