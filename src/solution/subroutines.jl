
function benchmark_subroutines(sol::Solution{T}, ode_wrap!::ODEWrapperState,
                               options::SolverOptions{T};
                               n_sample::Int64 = 1000) where T <: AbstractFloat

    @unpack nt, ny, np = get_dimensions(sol)
    @unpack method, stage_finder = options

    @assert 2 <= n_sample <= nt
    @assert nt > 1

    t_idxs = round.(Int64, LinRange(1, nt, n_sample))

    y = zeros(T, ny)
    f = zeros(T, ny)

    S2_runtime = 0.0
    FE_runtime = 0.0
    BC_runtime = 0.0

    @unpack stages = method

    sol_t = Vector{T}()
    sol_y = Vector{T}()
    sizehint!(sol_t, 1)
    sizehint!(sol_y, ny)
    # 32.4 ms - 28.8ms = 3.6 ms

    for n in t_idxs
        t = sol.t[n]
        y .= view(sol.y, 1+(n-1)*ny:n*ny)

        S2_stat = @timed begin
            append!(sol_t, t)
            append!(sol_y, y)
        end
        S2_runtime += S2_stat.time
        empty!(sol_t)
        empty!(sol_y)

        FE_stat = @timed ode_wrap!(f, t, y)
        FE_runtime += FE_stat.time

        # not this simple obviously
        BC_1 = @timed @.. f = y
        BC_2 = @timed @.. y = y + f
        BC_runtime += BC_1.time
        # BC_runtime += stages*BC_2.time
    end
    S2_runtime *= nt / length(t_idxs)
    FE_runtime *= sol.FE[1] / length(t_idxs)
    BC_runtime *= sol.FE[1] / length(t_idxs)

    println("")
    println("  Subroutine runtimes (seconds)  ")
    println("---------------------------------")
    println("save solution        | $(round(S2_runtime, sigdigits = 4))")
    println("function evaluations | $(round(FE_runtime, sigdigits = 4))")
    println("broadcast operations | $(round(BC_runtime, sigdigits = 4))")
    println("")
    return nothing
end