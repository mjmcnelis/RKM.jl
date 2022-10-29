
function get_precision(method::RungeKutta)
    # TODO: gives 1 allocation, do I need to preallocate matrix?
    typeof(method.butcher[1,1])
end
