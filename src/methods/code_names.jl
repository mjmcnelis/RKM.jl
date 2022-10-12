
function code_names(::RungeKutta, ::AdaptiveStepSize, ::Explicit)
    Dict(:Euler_1 => :E1,)
end
function code_names(::RungeKutta, ::Embedded, ::Explicit)
    Dict(:Fehlberg_12 => :F12,)
end

# nothing yet
function get_code_name(method)
    nothing
end