"""
    plot_ode(sol::Solution, method::ODEMethod, plot::Function)

Plots the ODE solution `y(t)` from `sol` with the function `plot` (usually `Plots.plot`).

Note: pass `plot` as a function to reduce RKM's precompilation time.
"""
function plot_ode(sol::Solution, method::ODEMethod, plot::Function)
    # TODO: put some more functionality (panel?)
    y, t = get_solution(sol)
    
    @unpack code_name = method 

    plot(t, y; size = (900, 600), linewidth = 2, legendtitle = code_name,
         legend = :outertopright, legendtitlefontsize = 12, legendfontsize = 12,
         ylabel = "y", yguidefontsize = 14, ytickfontsize = 12,
         xlabel = "t", xguidefontsize = 14, xtickfontsize = 12,) |> display

    # dt = t[2:end] .- t[1:end-1]
    # plot(t[1:end-1], dt; size = (900, 600), linewidth = 2,
    #      label = code_name, legend = :outertopright, legendfontsize = 12,
    #      ylabel = "Δt", yguidefontsize = 14, ytickfontsize = 12,
    #      xlabel = "t", xguidefontsize = 14, xtickfontsize = 12,) |> display
    nothing
end
