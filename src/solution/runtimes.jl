struct SolverRuntimes
    """Runtime of time evolution loop [seconds]"""
    evolution_time::MVector{1,Float64}
    """Function evaluation runtime [seconds]"""
    FE_time::MVector{1,Float64}
    """State Jacobian runtime [seconds]"""
    Jy_time::MVector{1,Float64}
    """Parameter Jacobian runtime [seconds]"""
    Jp_time::MVector{1,Float64}
    """Jacobian-vector product runtime [seconds]"""
    Jv_time::MVector{1,Float64}
    """Linear solve runtime [seconds]"""
    LS_time::MVector{1,Float64}
    """Save solution runtime [seconds]"""
    save_time::MVector{1,Float64}
end

function SolverRuntimes()
    evolution_time = MVector{1,Float64}(0.0)
    FE_time = MVector{1,Float64}(0.0)
    Jy_time = MVector{1,Float64}(0.0)
    Jp_time = MVector{1,Float64}(0.0)
    Jv_time = MVector{1,Float64}(0.0)
    LS_time = MVector{1,Float64}(0.0)
    save_time = MVector{1,Float64}(0.0)

    return SolverRuntimes(evolution_time, FE_time, Jy_time, Jp_time,
                          Jv_time, LS_time, save_time)
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
    runtimes.Jy_time .= 0.0
    runtimes.Jp_time .= 0.0
    runtimes.Jv_time .= 0.0
    runtimes.LS_time .= 0.0
    runtimes.save_time .= 0.0
    return nothing
end

function compute_runtimes!(runtimes::SolverRuntimes, config::RKMConfig,
                           loop_stats::NamedTuple, save_time::Vector{Float64})

    ode_wrap_y! = config.ode_wrap_y!
    # only used for parameter jacobian but should add to FE in general stats
    # ode_wrap_p! = config.ode_wrap_p!
    state_jacobian = config.state_jacobian
    sensitivity = config.sensitivity

    root_finder = config.root_finder

    evolution_time = runtimes.evolution_time
    FE_time = runtimes.FE_time
    Jy_time = runtimes.Jy_time
    Jp_time = runtimes.Jp_time
    Jv_time = runtimes.Jv_time
    LS_time = runtimes.LS_time

    evolution_time[1] = loop_stats.time
    FE_time[1] = ode_wrap_y!.runtime[1]

    Jy_time[1] = state_jacobian.runtime[1]

    if hasproperty(sensitivity, :param_jacobian)
        param_jacobian = sensitivity.param_jacobian
        Jp_time[1] = param_jacobian.runtime[1]
    end

    if hasproperty(sensitivity, :jacobian_vector)
        jacobian_vector = sensitivity.jacobian_vector
        Jv_time[1] = jacobian_vector.runtime[1]
    end

    if root_finder isa Newton
        LS_time[1] = root_finder.runtime[1]
    else
        LS_time[1] = 0.0
    end

    runtimes.save_time[1] = save_time[1]

    return nothing
end

function get_subroutine_times(runtimes::SolverRuntimes)
    if iszero(runtimes.FE_time[1])
        println("")
        @warn "No subroutine times available (set time_subroutine = true)"

        FE_time = "N/A"
        Jy_time = "N/A"
        Jp_time = "N/A"
        Jv_time = "N/A"
        LS_time = "N/A"
        save_time = "N/A"
    else
        FE_time = round(runtimes.FE_time[1], sigdigits = 4)
        Jy_time = round(runtimes.Jy_time[1], sigdigits = 4)
        Jp_time = round(runtimes.Jp_time[1], sigdigits = 4)
        Jv_time = round(runtimes.Jv_time[1], sigdigits = 4)
        LS_time = round(runtimes.LS_time[1], sigdigits = 4)
        save_time = round(runtimes.save_time[1], sigdigits = 4)
    end

    println("")
    println("Subroutine times (seconds)")
    println("---------------------------------")
    println("function evaluations | $FE_time")
    println("state jacobian       | $Jy_time")
    println("parameter jacobian   | $Jp_time")
    println("jacobian-vector      | $Jv_time")
    println("linear solve         | $LS_time")
    println("save solution        | $save_time")
    return nothing
end
