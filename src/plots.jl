
function plot_ode(sol, method, plot::Function)

    y, t = reshape_solution(sol)

    plot(t, y; size = (900, 600), linewidth = 2, legendtitle = method.code_name,
         legend = :outertopright, legendtitlefontsize = 12, legendfontsize = 12,
         ylabel = "y", yguidefontsize = 14, ytickfontsize = 12,
         xlabel = "t", xguidefontsize = 14, xtickfontsize = 12,) |> display
#     plot(t[1:end-1], dt; size = (900, 600), linewidth = 2,
#          label = method.code_name, legend = :outertopright, legendfontsize = 12,
#          ylabel = "Δt", yguidefontsize = 14, ytickfontsize = 12,
#          xlabel = "t", xguidefontsize = 14, xtickfontsize = 12,) |> display
end
