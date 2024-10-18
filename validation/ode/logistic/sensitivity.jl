using Revise, RKM
using OrdinaryDiffEq, SciMLSensitivity
using ForwardDiff: Dual, value, partials
using LinearSolve, Setfield
using Plots; plotly()
!(@isdefined dy_dt!) ? include("$RKM_root/validation/ode/logistic/equations.jl") : nothing
include("$RKM_root/validation/ode/logistic/parameters.jl")

options = SolverOptions(options)

t0 = -5.0
tf = 5.0
dt0 = 0.1

N = 2
p = [0.5 - 0.25*(i-1.0)/(N-1.0+eps(1.0)) for i in 1:N]
y0 = Float64[]
for i = eachindex(p)
    push!(y0, exp(t0) / (1.0 + exp(t0)) - p[i])
end

@time sol = evolve_ode(y0, t0, tf, dt0, dy_dt!, options, p)
S, t = get_sensitivity(sol)

ny = sol.dimensions[1]
nt = length(sol.t)
np = length(p)

# time slice: dy_1/dp1 dy2/dp1 dy1/dp2 dy2/dp2

# get time slice (transpose?)
# @time a = reshape(view(S, 2, :), ny, np)

# get parameter slice
# @time b = view(S, :, 1:ny)

# get_stats(sol)

GC.gc()
println("")

#------------------------------------------
code_name = options.method.code_name
if code_name == "BE1"
    alg = ImplicitEuler(; autodiff = false)
elseif code_name == "TRBDF2"
    alg = TRBDF2(; autodiff = false)
elseif code_name == "RK4"
    alg = RK4()
end
tspan = (t0, tf)

# Dual numbers
zers = zeros(N)
p_dual = Vector{Dual{Nothing, Float64, N}}()
for i in 1:N
    zers[i] = 1.0
    push!(p_dual, Dual{Nothing}(p[i], zers...))
    zers[i] = 0.0
end
prob_dual = ODEProblem{true, SciMLBase.FullSpecialize}(f_ord, convert.(eltype.(p_dual), y0), tspan, p_dual)
@time sol_dual = solve(prob_dual, alg; dt = dt0, saveat = dt0, adaptive = false)
yv_dual = hcat([value.(sol_dual.u[i]) for i in eachindex(sol_dual.u)]...)'
yp_dual = hcat([vcat(partials.(sol_dual.u[i])...) for i in eachindex(sol_dual.u)]...)'

GC.gc()

# Forward Sensitivity (direct method)
prob_fs = ODEForwardSensitivityProblem(f_ord, y0, tspan, p)
@time sol_fs = solve(prob_fs, alg; dt = dt0, saveat = dt0, adaptive = false)
yv_fs = hcat(sol_fs.u...)'[:,1:N]
yp_fs = hcat(sol_fs.u...)'[:, N+1:end]

# peak is t = 0
t_idx = 51
a = reshape(view(S, t_idx, :), ny, np)
a_dual = reshape(view(yp_dual, t_idx, :), ny, np)
a_fs = reshape(view(yp_fs, t_idx, :), ny, np)

if N <= 2
    plot(t[1:end-1], S[1:end-1,:]) |> display
    plot(sol_dual.t, yp_dual) |> display
    plot(sol_fs.t, yp_fs) |> display
end

println("\ndone")