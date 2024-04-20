
# note: thread on DataType kwarg type-instability
# https://discourse.julialang.org/t/type-stability-when-datatype-is-passed-as-a-keyword-argument/52770

# TODO: not sure if something like this is needed (try again later)
# abstract type RKMSolution end
# struct Solution{T} <: RKMSolution where T <: AbstractFloat

"""
Stores the solution vector `y(t)` of the ODE system in linear column format.

For example, the solution set `{y(0.0) = [1.0, 2.0, 3.0], y(0.5) = [4.0, 5.0, 6.0]}`
is stored as `y = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0], t = [0.0, 0.5]`.
"""
struct Solution{T <: AbstractFloat}
    """State vector of solution (stored as linear column)"""
    y::Vector{T}
    """Time vector of solution"""
    t::Vector{T}
    """Number of function evaluations"""
    FE::MVector{1,Int64}
    """Number of Jacobian evaluations"""
    JE::MVector{1,Int64}
    """Step rejection rate (percentage)"""
    rejection_rate::MVector{1,Float64}
    """Runtime of ODE solver (excludes configuration) [seconds]"""
    runtime::MVector{1,Float64}
    """Memory used to store solution set (t,y) [bytes]"""
    solution_size::MVector{1,Int64}
    """Memory used to configure the solver [bytes]"""
    config_memory::MVector{1,Int64}
    """Excess memory used in time evolution loop [bytes]"""
    excess_memory::MVector{1,Int64}
    """Number of dynamical variables (assumed to be fixed)"""
    dimensions::MVector{1,Int64}
    """Float precision type"""
    precision::Type{T}
end

"""
    Solution(precision::Type{T} = Float64) where T <: AbstractFloat

Outer constructor for `Solution`.

Required parameters: `precision`
"""
function Solution(precision::Type{T} = Float64) where T <: AbstractFloat
    y = Vector{precision}()
    t = Vector{precision}()
    FE = MVector{1,Int64}(0)
    JE = MVector{1,Int64}(0)
    rejection_rate = MVector{1,Float64}(0.0)
    runtime = MVector{1,Float64}(0.0)
    solution_size = MVector{1,Int64}(0)
    config_memory = MVector{1,Int64}(0)
    excess_memory = MVector{1,Int64}(0)
    dimensions = MVector{1,Int64}(0.0)

    return Solution(y, t, FE, JE, rejection_rate, runtime, solution_size,
                    config_memory, excess_memory, dimensions, precision)
end

function clear_solution!(sol::Solution)
    @unpack y, t, FE, JE, rejection_rate, runtime,
            solution_size, config_memory, excess_memory = sol
    empty!(y)
    empty!(t)
    sizehint!(y, 0)
    sizehint!(t, 0)
    FE .= 0
    JE .= 0
    rejection_rate .= 0.0
    runtime .= 0.0
    solution_size .= 0
    config_memory .= 0
    excess_memory .= 0
    return nothing
end

"""
    sizehint_solution!(sol::Solution, t0::T, tf::T, dt0::T,
                       dimensions::Int64) where T <: AbstractFloat

Applies `sizehint!` to the vector fields `y` and `t` in the solution `sol`.

Required parameters: `sol`, `t0`, `tf`, `dt0`, `dimensions`
"""
function sizehint_solution!(sol::Solution, t0::T, tf::T, dt0::T,
                            dimensions::Int64) where T <: AbstractFloat
    # TODO: whether or not I call this depends if save at regular intervals
    steps = round((tf - t0)/dt0) |> Int64
    # TODO: why is steps+1 not sufficient?
    sizehint!(sol.y, dimensions*(steps + 2))
    sizehint!(sol.t, steps + 2)
    return nothing
end

"""
    get_solution(sol::Solution)

Returns the solution tuple `(y,t)` from `sol`. The solution vector `y`, which has a length
`D*N`, is reshaped into an `N x D` matrix (`N` = time steps, `D` = dimensions).

For example, the solution set `{y(0.0) = [1.0, 2.0, 3.0], y(0.5) = [4.0, 5.0, 6.0]}`
is stored in linear column format as `y = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]`. The solution
vector is then reshaped as `y = [1.0 2.0 3.0; 4.0 5.0 6.0]`.
"""
function get_solution(sol::Solution)
    @unpack y, t, dimensions = sol
    # TODO: replace length(t) if use deleteat for PDEs
    y = reshape(y, dimensions[1], length(t)) |> transpose
    return y, t
end

function compute_stats!(sol::Solution, save_solution::Bool, adaptive::AdaptiveStepSize,
                        timer::TimeLimit, stage_finder::StageFinder,
                        loop_stats::NamedTuple, config_bytes::Int64)

    @unpack jacobian_method = stage_finder
    @unpack rejection_rate, JE, runtime, solution_size, config_memory, excess_memory = sol

    rejection_rate .= compute_step_rejection_rate(adaptive, timer)
    JE .= jacobian_method.evaluations
    runtime .= loop_stats.time
    solution_size .= sizeof(sol.y) + sizeof(sol.t)
    config_memory .= config_bytes
    excess_memory .= loop_stats.bytes
    if !(save_solution && adaptive isa Fixed)
        excess_memory .-= solution_size
    end
    return nothing
end

function get_stats(sol::Solution)
    @unpack t, rejection_rate, FE, JE, runtime,
            solution_size, config_memory, excess_memory = sol
    println("time steps           = $(length(t))")
    println("step rejection rate  = $(round(rejection_rate[1], sigdigits = 4)) %")
    println("function evaluations = $(FE[1])")
    println("jacobian evaluations = $(JE[1])")
    println("evolution runtime    = $(round(runtime[1], sigdigits = 4)) seconds")
    println("solution size        = $(format_bytes(solution_size[1]))")
    println("config memory        = $(format_bytes(config_memory[1]))")
    println("excess memory        = $(format_bytes(excess_memory[1]))")
end
