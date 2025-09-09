struct SolverRuntimes
    """Runtime of time evolution loop [seconds]"""
    evolution_time::MVector{1,Float64}
    """Function evaluation runtime [seconds]"""
    FE_time::MVector{1,Float64}
    """Jacobian evaluation runtime [seconds]"""
    JE_time::MVector{1,Float64}
    """Linear solve runtime [seconds]"""
    LS_time::MVector{1,Float64}
    """Save solution runtime [seconds]"""
    save_time::MVector{1,Float64}
end

function SolverRuntimes()
    evolution_time = MVector{1,Float64}(0.0)
    FE_time = MVector{1,Float64}(0.0)
    JE_time = MVector{1,Float64}(0.0)
    LS_time = MVector{1,Float64}(0.0)
    save_time = MVector{1,Float64}(0.0)

    return SolverRuntimes(evolution_time, FE_time, JE_time, LS_time, save_time)
end

function Base.show(io::IO, runtimes::SolverRuntimes)
    println("")
    evolution_time = runtimes.evolution_time
    println("Evolution runtime = $(round(evolution_time[1], sigdigits = 4)) seconds")
    get_subroutine_times(runtimes)

    return nothing
end

function clear_runtimes!(runtimes::SolverRuntimes)
    runtimes.evolution_time .= 0.0
    runtimes.FE_time .= 0.0
    runtimes.JE_time .= 0.0
    runtimes.LS_time .= 0.0
    runtimes.save_time .= 0.0
    return nothing
end

function compute_runtimes!(runtimes::SolverRuntimes, config::RKMConfig,
                           loop_stats::NamedTuple, save_time::Vector{Float64})

    ode_wrap_y! = config.ode_wrap_y!
    state_jacobian = config.state_jacobian
    root_finder = config.root_finder

    evolution_time = runtimes.evolution_time
    FE_time = runtimes.FE_time
    JE_time = runtimes.JE_time
    LS_time = runtimes.LS_time

    evolution_time[1] = loop_stats.time
    FE_time[1] = ode_wrap_y!.runtime[1]
    JE_time[1] = state_jacobian.runtime[1]
    if root_finder isa Newton
        LS_time[1] = root_finder.runtime[1]
    else
        LS_time[1] = 0.0
    end
    runtimes.save_time[1] = save_time[1]

    return nothing
end

function get_subroutine_times(runtimes::SolverRuntimes)
    println("")
    println("Subroutine times (seconds)")
    println("---------------------------------")
    println("function evaluations | $(round(runtimes.FE_time[1], sigdigits = 4))")
    println("jacobian evaluations | $(round(runtimes.JE_time[1], sigdigits = 4))")
    println("linear solve         | $(round(runtimes.LS_time[1], sigdigits = 4))")
    println("save solution        | $(round(runtimes.save_time[1], sigdigits = 4))")
    return nothing
end
