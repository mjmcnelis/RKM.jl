# TODO: move somewhere else
struct JacobianException
    msg::String
end
Base.showerror(io::IO, e::JacobianException) = print(io, "JacobianException: ", e.msg)
function jacobian_error(args...; kwargs...) 
    msg = "using implicit method but no jacobian has been specified or computed"
    throw(JacobianException(msg))
end

function evolve_ode(y0, dy_dt!::Function; jacobian!::Function = jacobian_error, # TEMP
                                          parameters::Parameters)

    @unpack adaptive, method, t_span, timer = parameters
    @unpack t0, tf, dt0 = t_span
    
    @unpack stages, precision, iteration = method 

    # initial conditions
    y  = precision[copy(y0)...]
    t  = MVector{1}(t0)
    dt = MVector{2}(dt0, dt0)

    dimensions = size(y, 1)
    
    # note: should not be SA in general but still may want option if size small
    # note: keep in mind of ForwardDiff issues we had with PaT
    dy    = zeros(precision, stages, dimensions) 
    y_tmp = zeros(precision, dimensions)
    f_tmp = zeros(precision, dimensions)
    f     = zeros(precision, dimensions)
    # f = @SVector zeros(precision, dimension)  # want to test it out though
   
    # TEMP for step doubling (embedded too probably)
    y1 = zeros(precision, dimensions)
    y2 = zeros(precision, dimensions)
    error = zeros(precision, dimensions)

    # initalize solution
    sol = Solution(; precision, dimensions) 

    while true
        update_solution!(sol, y, t)
      
        # TODO: see if can pass kwargs
        evolve_one_time_step!(method, iteration, adaptive, y, t, dt, dy_dt!, 
                              dy, y_tmp, f_tmp, f, y1, y2, error, jacobian!)

        check_time(t, tf, timer) || break
        t .+= dt[1]
    end
    sol
end

function update_solution!(sol, y::Vector{T}, t::MVector{1,T}) where T <: AbstractFloat
    # push!(sol.y, copy(y))
    for i in eachindex(y) 
        append!(sol.y[i], y[i])
    end
    append!(sol.t, t)
    nothing 
end
