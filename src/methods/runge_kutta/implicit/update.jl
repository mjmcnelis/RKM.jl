
# TODO: so far, routine only works for an explicit, primary method
function fixed_runge_kutta_step!(method::RungeKutta, iteration::DiagonalImplicit, 
                                 y::Vector{<:AbstractFloat}, 
                                 t::Float64, dt::Float64, dy_dt!::Function, 
                                 dy::Matrix{<:AbstractFloat}, 
                                 y_tmp::Vector{<:AbstractFloat}, 
                                 f_tmp::Vector{<:AbstractFloat})
           
    butcher = method.butcher                          
    ncol   = size(butcher, 2)  
    stages = ncol - 1       
    
    root_solver = "fixed_point"        # will just use fixed point iteration for now 
    eps_root = 1e-8 
    max_iterations = 10

    c = view(butcher, 1:ncol-1, 1)                   
    A = view(butcher, 1:ncol-1, 2:ncol)             
    b = view(butcher, ncol, 2:ncol)      

    for i = 1:stages  
        t_tmp = t + c[i]*dt 

        # sum over known stages 
        y_tmp .= y 
        for j = 1:i-1 
            for k in eachindex(y_tmp)
                y_tmp[k] += A[i,j] * dy[j,k]
            end
        end 

        # will just do fixed point iteration w/o any breaks for now 
        for n = 1:max_iterations
            for k in eachindex(y_tmp)
                y_tmp[k] += A[i,i] * dy[i,k]
            end 
            dy_dt!(f_tmp, t_tmp, y_tmp)
            # TEMP undo
            for k in eachindex(y_tmp)
                y_tmp[k] -= A[i,i] * dy[i,k]
            end
            dy[i,:] .= dt .* f_tmp 
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

function evolve_one_time_step!(method::RungeKutta, iteration::DiagonalImplicit,
                               adaptive::Fixed,
                               y::Vector{<:AbstractFloat}, t::Vector{Float64}, 
                               dt::Vector{Float64}, dy_dt!::Function, 
                               dy::Matrix{<:AbstractFloat}, y_tmp::Vector{<:AbstractFloat}, 
                               f_tmp::Vector{<:AbstractFloat}, f::Vector{<:AbstractFloat},
                               args...) 
    # TODO: not sure why putting dy_dt! here this kills allocations
    # costs an extra stage but saves on allocations
    dy_dt!(f, t[1], y)

    fixed_runge_kutta_step!(method, iteration, y, t[1], dt[1], dy_dt!, dy, y_tmp, f_tmp)
    y .= y_tmp
    nothing
end
