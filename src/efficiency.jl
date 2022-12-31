
function efficiency_curve(y0::Union{T, Vector{T}}, y_exact::Function, dy_dt!::Function; 
             precision::Type{T2}, methods::OrderedDict{<:AdaptiveStepSize, <:Vector}, 
             epsilon_vect::Vector{Float64}, t_range::TimeRange, plot::Function, 
             plot!::Function) where {T <: AbstractFloat, T2 <: AbstractFloat}
    plt = plot()
    for key in keys(methods)
            adaptive = key
        for method in methods[key]
            FE = Vector{Float64}()
            mean_err = Vector{Float64}()

            # TODO: sort out how to do efficiency for fixed time step
            for epsilon in epsilon_vect
                # @show epsilon
                adaptive = @set adaptive.epsilon = epsilon 
                parameters = Parameters(; adaptive, method, t_range)

                sol = evolve_ode(y0, dy_dt!; parameters, precision)
                y, t = get_solution(sol)

                y_ex = zeros(Double64, size(y)...)
                err = zeros(Double64, length(t))

                # TODO: how to skip passing dimensions?
                for i in eachindex(t)
                    y_ex[i,:] = y_exact(t[i], sol.dimensions)
                    err[i] = norm(y[i,:] .- y_ex[i,:], adaptive.p_norm)
                end
                
                append!(FE, sol.FE)
                append!(mean_err, mean(err))
            end
            code_name = method.code_name*adaptive_code_label(adaptive)
            linestyle = get_linestyle(adaptive)

            # TODO: separate legends for SDRK, ERK
            plot!(FE, mean_err; 
                  size = (900, 600), linewidth = 2, linestyle, label = code_name,
                  legend = :outertopright, legendtitlefontsize = 12, legendfontsize = 12,
                  ylabel = "Mean norm error", yguidefontsize = 14, ytickfontsize = 12,
                  xlabel = "Function evaluations", xguidefontsize = 14, xtickfontsize = 12,
                #   ylims = (1e-30, 1e0), xlims = (1e2, 1e7), 
                  ylims = (1e-15, 1e0), xlims = (1e2, 1e5), 
                  xaxis = :log, yaxis = :log)
        end
    end
    plt
end