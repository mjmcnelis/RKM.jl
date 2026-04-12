using Plots; plotly()
using LinearAlgebra: I, tr

function x_exact(t)
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

function midpoint_rule(x0, t_vect; discrete = true)
    nx = length(x0)
    nt = length(t_vect)

    t0 = t_vect[1]
    t_idxs = 2:length(t_vect)
    t_prev = t0

    x = copy(x0)
    x_tmp = copy(x0)
    xM = [x0...]

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
        append!(xM, x)
    end
    xM = reshape(xM, nx, nt) |> transpose
    return xM
end

function scalar_generator(x0, t_vect; discrete = true, order = 1)
    nx = length(x0)
    nt = length(t_vect)

    t0 = t_vect[1]
    t_idxs = 2:length(t_vect)
    t_prev = t0

    x = copy(x0)
    x_tmp = copy(x0)
    xG = [x0...]

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
        append!(xG, x)
    end
    xG = reshape(xG, nx, nt) |> transpose
    return xG
end

function matrix_generator(x0, t_vect; discrete = true, order = 1)
    nx = length(x0)
    nt = length(t_vect)

    t0 = t_vect[1]
    t_idxs = 2:length(t_vect)
    t_prev = t0

    x = copy(x0)
    x_tmp = copy(x0)
    xG = [x0...]

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
        append!(xG, x)
    end
    xG = reshape(xG, nx, nt) |> transpose
    return xG
end

dt = 0.5
t0 = 0.0
tf = 5.0
x0 = [0.0, 1.0]
t_vect = t0:dt:tf

# x_SG = scalar_generator(x0, t_vect; order = 1)
x_MG = matrix_generator(x0, t_vect; order = 1)
x_MR = midpoint_rule(x0, t_vect)
# x = x_SG
x = x_MG
# x = x_MR

t_fine = t0:0.1:tf
plt = plot(t_fine, x_exact(t_fine), color = [:black :gray], linewidth = 1.5);
plot!(t_vect, x, color = [:red :blue], linewidth = 1.5, line = :dash)
scatter!(t_vect, x, color = [:red :blue], markersize = [4,4])
display(plt)

sum(abs2, x_SG .- x_exact(t_vect)) |> display
sum(abs2, x_MG .- x_exact(t_vect)) |> display
sum(abs2, x_MR .- x_exact(t_vect)) |> display

# plt = plot(t_vect, x_SG .- x_exact(t_vect), color = [:red :blue],
#            label = ["Δx_1 (sg)";; "Δx_2 (sg)"], linewidth = 1.5)
# plot!(t_vect, x_MG .- x_exact(t_vect), color = [:red :blue],
#       label = ["Δx_1 (mg)";; "Δx_2 (mg)"], linewidth = 1.5, line = :dash)
# plot!(t_vect, x_MR .- x_exact(t_vect), color = [:red :blue],
#       label = ["Δx_1 (mid)";; "Δx_2 (mid)"], linewidth = 1.5, line = :dot)
# display(plt)

println("\ndone")