
# TODO: so far, routine only works for an explicit, primary method
# TODO: try using the where {T} notation 
function fixed_runge_kutta_step!(method::RungeKutta, iteration::DiagonalImplicit, 
                                 y::Vector{<:AbstractFloat}, 
                                 t::Float64, dt::Float64, dy_dt!::Function, 
                                 dy::Matrix{<:AbstractFloat}, 
                                 y_tmp::Vector{<:AbstractFloat}, 
                                 f_tmp::Vector{<:AbstractFloat}, jacobian!::Function)
           
    butcher = method.butcher                          
    ncol   = size(butcher, 2)  
    stages = ncol - 1       
    
    # root_solver = "fixed_point"        # will just use fixed point iteration for now 
    root_solver = "newton_fast"
    eps_root = 1e-8 
    max_iterations = 10

    c = view(butcher, 1:ncol-1, 1)                   
    A = view(butcher, 1:ncol-1, 2:ncol)             
    b = view(butcher, ncol, 2:ncol)  

    # TEMP 
    L = length(y) 
    J = zeros(L, L)             # allocates

    for i = 1:stages  
        # evaluate jacobian
        if root_solver == "newton_fast" # newton fast 
            jacobian!(J, t, y)        
            J .*= (-A[i,i]*dt)
            for i in eachindex(y) 
                J[i,i] += 1.0
            end
        end

        t_tmp = t + c[i]*dt 
        # sum over known stages 
        y_tmp .= y 
        for j = 1:i-1 
            for k in eachindex(y_tmp)
                y_tmp[k] += A[i,j] * dy[j,k]
            end
        end 

        # TEMP iterate w/o any breaks for now 
        for n = 1:max_iterations
            for k in eachindex(y_tmp)
                y_tmp[k] += A[i,i] * dy[i,k]
            end 
            dy_dt!(f_tmp, t_tmp, y_tmp)
            # TEMP undo addition (for minimizing allocations)
            for k in eachindex(y_tmp)
                y_tmp[k] -= A[i,i] * dy[i,k]
            end

            if root_solver == "fixed_point"
                dy[i,:] .= dt .* f_tmp 
            elseif root_solver == "newton_fast"
                # TODO: try to solve for dy directly instead of d(dy)
                for k in eachindex(f_tmp) 
                    f_tmp[k] = dy[i,k] - dt*f_tmp[k]
                end
                # from python 
                # g = z - dt*y_prime(t + dt*c[i], y + dy + z*Aii)
                dy[i, :] .-= J \ f_tmp
            end
        end
    end    

    y_tmp .= y                                    
    for j = 1:stages 
        for i in eachindex(y_tmp)
            y_tmp[i] += b[j] * dy[j,i]
        end
    end
    nothing
end

# TODO: use alternative indents
function evolve_one_time_step!(method::RungeKutta, iteration::DiagonalImplicit,
                               adaptive::Fixed,
                               y::Vector{<:AbstractFloat}, t::Vector{Float64}, 
                               dt::Vector{Float64}, dy_dt!::Function, 
                               dy::Matrix{<:AbstractFloat}, y_tmp::Vector{<:AbstractFloat}, 
                               f_tmp::Vector{<:AbstractFloat}, f::Vector{<:AbstractFloat},
                               y1, y2, error, jacobian!, args...) 
    # TODO: not sure why putting dy_dt! here this kills allocations
    # costs an extra stage but saves on allocations
    # TODO: want to ultimately remove this 
    dy_dt!(f, t[1], y)

    fixed_runge_kutta_step!(method, iteration, y, t[1], dt[1], dy_dt!, dy, y_tmp, f_tmp, 
                            jacobian!)
    y .= y_tmp
    nothing
end
