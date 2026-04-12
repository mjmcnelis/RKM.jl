using Plots; plotly()
using LinearAlgebra: I, tr

function f_calc(t)
    f = [0.5*t*exp(-t), exp(-t)]
    return f
end

function dfdt_calc(t)
    dfdt = [0.5*(1.0 - t)*exp(-t), -exp(-t)]
    return dfdt
end

function dfdt2_calc(t)
    dfdt2 = [0.5*(t - 2.0)*exp(-t), exp(-t)]
    return dfdt2
end

function A_calc(t)
    A = [-2.0 t
         0.0  -1.0]
    return A
end

function Ainv_calc(t)
    Ainv = [-0.5 -0.5*t
            0.0  -1.0]
    return Ainv
end

function dAinvdt_calc(t)
    dAinvdt = [0.0 -0.5
               0.0  0.0]
    return dAinvdt
end

function z0_calc(t, t_prev)
    z0 = [2.0*(t-t_prev) 0.5*(t_prev^2-t^2)
          0.0            t-t_prev]
    return z0
end

function x_exact(t)
    x = zeros(length(t), 2)
    x[:,1] .= (t.^2 .- t .+ 1.0).*exp.(-t) - exp.(-2.0.*t)
    x[:,2] = (t .+ 1.0).*exp.(-t)
    return x
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

        Ainv = Ainv_calc(t)
        dAinvdt = dAinvdt_calc(t)
        f = f_calc(t)
        dfdt = dfdt_calc(t)
        dfdt2 = dfdt2_calc(t)
        z0 = z0_calc(t, t_prev)

        xG0 = f + exp(-z0)*(x_tmp - f)
        xG1 = (I - exp(-z0)*(I + z0))*Ainv*dfdt
        # note: this does help even though it's ignoring commutators
        xG2 = (I - exp(-z0)*(I + z0 + 0.5*z0^2))*Ainv*(dAinvdt*dfdt + Ainv*dfdt2)
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

xG = scalar_generator(x0, t_vect; discrete = true, order = 1)
# xG = matrix_generator(x0, t_vect; discrete = true, order = 1)

t_fine = t0:0.1:tf
plt = plot(t_fine, x_exact(t_fine), color = [:black :gray], linewidth = 1.5);
plot!(t_vect, xG, color = [:red :blue], linewidth = 1.5, line = :dash)
scatter!(t_vect, xG, color = [:red :blue], markersize = [4,4])

display(plt)

println("\ndone")