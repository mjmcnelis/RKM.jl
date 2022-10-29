
using Revise
using RKM

# TODO: better interface? (MTK?)
function y_prime(t, y) 
    (y + 0.5) * (0.5 - y)
end

# should handle scalar and vector cases? or just vector
y0 = -1.0
t_span = TimeSpan(; t0 = -10.0, tf = 10.0, dt0 = 2.0)

adaptive   = Fixed()
method     = Heun2()
parameters = Parameters(; adaptive, method)


@time sol = evolve_ode(y0, t_span, y_prime; parameters)

# @show sol.t sol.dt

println("\ndone")

