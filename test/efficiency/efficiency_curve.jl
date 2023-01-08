import RKM: AdaptiveStepSize, adaptive_code_label, get_linestyle
import LinearAlgebra: norm
import Statistics: mean
import Plots: plot, plot!
import Setfield: @set

function efficiency_curve(y0::Union{T, Vector{T}}, y_exact::Function, dy_dt!::Function; 
             precision::Type{T2}, methods::OrderedDict{<:AdaptiveStepSize, <:Vector}, 
             epsilon_vect::Vector{Float64}, t_range::TimeRange, xlims = (1e2, 1e5),
             ylims = (1e-20, 1e0)) where {T <: AbstractFloat, T2 <: AbstractFloat}
    plt = plot()

    for key in keys(methods)
            adaptive = key
        for method in methods[key]
            FE = Vector{Float64}()
            mean_err = Vector{Float64}()
            p_norm = adaptive.p_norm

            # TODO: sort out how to do efficiency for fixed time step
            for epsilon in epsilon_vect
                adaptive = @set adaptive.epsilon = epsilon 
                parameters = Parameters(; adaptive, method, t_range)

                sol = evolve_ode(y0, dy_dt!; parameters, precision)
                y, t = get_solution(sol)

                # TODO: redo this for better memory efficiency
                y_ex = zeros(BigFloat, size(y)...)
                err = zeros(BigFloat, length(t))

                # TODO: how to skip passing dimensions?
                for i in eachindex(t)
                    i == 1 ? continue : nothing # skip initial condition
                    y_ex[i,:] = y_exact(t[i], sol.dimensions)
                    err[i] = norm(y[i,:] .- y_ex[i,:], p_norm)
                end
                
                append!(FE, sol.FE)
                append!(mean_err, mean(err))

                sol.FE[1] < xlims[2] || break
                if !isnan(mean(err))
                    mean(err) > ylims[1] || break
                end
            end
            code_name = method.code_name*adaptive_code_label(adaptive)
            linestyle = get_linestyle(adaptive)
            norm_label = isinf(p_norm) ? "∞" : "$(Int(p_norm))"

            # TODO: separate legends for SDRK, ERK
            plot!(FE, mean_err; 
                  size = (900, 600), linewidth = 2, linestyle, label = code_name,
                  legend = :outertopright, legendtitlefontsize = 12, legendfontsize = 12,
                  xlabel = "Function evaluations", ylabel = "Mean L$norm_label error",
                  yguidefontsize = 14, ytickfontsize = 12, yaxis = :log, ylims,
                  xguidefontsize = 14, xtickfontsize = 12, xaxis = :log, xlims,
                )
        end
    end
    plt
end