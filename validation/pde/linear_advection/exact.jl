
    # TODO: construct exact numerical solution in another file
    # try LxF 
    # a = 1.0
    # dx = 0.1 
    # C = a*dt[1]/dx
    # L = length(y) 
    # A = zeros(L, L)
#=
    # backward central (Neumann BC)
    A[1,1] = 1.0 - C/2.0
    A[1,2] = C/2.0
    A[end,end-1] = -C/2.0
    A[end,end] = 1.0 + C/2.0
    for i = 2:L-1 
        A[i,i-1] = -C/2.0 
        A[i,i]   = 1.0
        A[i,i+1] = C/2.0
    end
=#

    # backup backward LxF (Neumann BC)
#=
    A[1,1] = (3.0 - C)/2.0
    A[1,2] = (C - 1.0)/2.0
    A[end,end-1] = -(C + 1.0)/2.0
    A[end,end] = (C + 3.0)/2.0
    for i = 2:L-1 
        A[i,i-1] = -(C + 1.0)/2.0 
        A[i,i]   = 2.0
        A[i,i+1] = (C - 1.0)/2.0
    end
=#
    # y .= A \ y