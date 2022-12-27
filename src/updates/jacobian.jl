
struct JacobianException
    msg::String
end

Base.showerror(io::IO, e::JacobianException) = print(io, "JacobianException: ", e.msg)

function jacobian_error(args...; kwargs...)
    msg = "using implicit method but no jacobian has been specified or computed"
    throw(JacobianException(msg))
end