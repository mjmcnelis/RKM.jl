# TODO: move somewhere else
struct JacobianException
    msg::String
end
Base.showerror(io::IO, e::JacobianException) = print(io, "JacobianException: ", e.msg)
function jacobian_error(args...; kwargs...)
    msg = "using implicit method but no jacobian has been specified or computed"
    throw(JacobianException(msg))
end

function evolve_ode(y0::Union{T, Vector{T}}, dy_dt!::F;
                    jacobian!::J, parameters::P) where {T <: AbstractFloat, F <: Function,
                                                        J <: Function, P <: Parameters}

    @unpack adaptive, method, t_span, timer = parameters
    @unpack t0, tf, dt0 = t_span

    @unpack stages, precision, iteration = method

    dimensions = size(y0, 1)

    # initial conditions
    y  = y0 .|> precision
    y isa Vector ? nothing : y = [y]

    t  = MVector{1,Float64}(t0)
    dt = MVector{2,Float64}(dt0, dt0)

    # note: should not be SA in general but still may want option if size small
    # note: keep in mind of ForwardDiff issues we had with PaT
    dy    = zeros(precision, dimensions, stages)
    y_tmp = zeros(precision, dimensions)
    f_tmp = zeros(precision, dimensions)
    f     = zeros(precision, dimensions)

    # TEMP for step doubling (embedded too probably)
    y1 = zeros(precision, dimensions)
    y2 = zeros(precision, dimensions)
    error = zeros(precision, dimensions)

    # initalize solution
    sol = Solution(; precision, dimensions)
    @unpack FE = sol

    sizehint_solution!(sol, t_span, dimensions)

    while true
        update_solution!(sol, y, t)

        # TODO: see if can pass kwargs
        evolve_one_time_step!(method, iteration, adaptive, FE, y, t, dt, dy_dt!,
                              dy, y_tmp, f_tmp, f, y1, y2, error, jacobian!)

        check_time(t, tf, timer) || break
        @.. t += dt[1]
    end
    sol
end
