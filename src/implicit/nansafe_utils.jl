#=
# requires Combinatorics package
function max_nan(args...)
    return sum(p -> maximum(p), permutations(args)) / factorial(length(args))
end
function min_nan(args...)
    return sum(p -> minimum(p), permutations(args)) / factorial(length(args))
end
=#

# TODO: allocates if mix float and dual numbers
function max_nan(args...)
    # idx = findfirst(x -> isnan(x), args)
    # if !isnothing(idx)
    #     ans = zero(args[idx])
    if any(x -> isnan(x), args)
        ans = zero(args[1])
        n = 0
        for i in eachindex(args)
            if isnan(args[i])
                ans += args[i]
                n += 1
            end
        end
        return ans/n
    else
        return max(args...)
    end
end

function min_nan(args...)
    if any(x -> isnan(x), args)
        ans = zero(args[1])
        n = 0
        for i in eachindex(args)
            if isnan(args[i])
                ans += args[i]
                n += 1
            end
        end
        return ans/n
    else
        return min(args...)
    end
end

function test_nansafe(; x = Dual(NaN, 1.0), y = Dual(NaN, 0.0),
                        b = 5.0, n = 1, r = RoundNearest)
    xv, xp = x.value, x.partials[1]
    yv, yp = y.value, y.partials[1]
    @assert isnan(xv)
    @assert iszero(xp) || isone(xp)
    @assert isnan(yv)
    @assert iszero(yp) || isone(yp)

    @show x
    println("")
    println("1.0 + x     = $(1.0 + x)")
    println("1.0 - x     = $(1.0 - x)")
    println("sqrt(x)     = $(sqrt(x))")
    println("cbrt(x)     = $(cbrt(x))")
    println("abs2(x)     = $(abs2(x))")
    println("inv(x)      = $(inv(x))")
    println("log(x)      = $(log(x))")
    println("log10(x)    = $(log10(x))")
    println("log2(x)     = $(log2(x))")
    println("log1p(x)    = $(log1p(x))")
    println("exp(x)      = $(exp(x))")
    println("exp2(x)     = $(exp2(x))")
    println("exp10(x)    = $(exp10(x))")
    println("expm1(x)    = $(expm1(x))")
    println("sin(x)      = $(sin(x))")
    println("cos(x)      = $(cos(x))")
    println("tan(x)      = $(tan(x))")
    println("sec(x)      = $(sec(x))")
    println("csc(x)      = $(csc(x))")
    println("cot(x)      = $(cot(x))")
    println("sind(x)     = $(sind(x))")
    println("cosd(x)     = $(cosd(x))")
    println("tand(x)     = $(tand(x))")
    println("secd(x)     = $(secd(x))")
    println("cscd(x)     = $(cscd(x))")
    println("cotd(x)     = $(cotd(x))")
    println("sinpi(x)    = $(sinpi(x))")
    println("cospi(x)    = $(cospi(x))")
    println("asin(x)     = $(asin(x))")
    println("acos(x)     = $(acos(x))")
    println("atan(x)     = $(atan(x))")
    println("asec(x)     = $(asec(x))")
    println("acsc(x)     = $(acsc(x))")
    println("acot(x)     = $(acot(x))")
    println("asind(x)    = $(asind(x))")
    println("acosd(x)    = $(acosd(x))")
    println("atand(x)    = $(atand(x))")
    println("asecd(x)    = $(asecd(x))")
    println("acscd(x)    = $(acscd(x))")
    println("acotd(x)    = $(acotd(x))")
    println("sinh(x)     = $(sinh(x))")
    println("cosh(x)     = $(cosh(x))")
    println("tanh(x)     = $(tanh(x))")
    println("sech(x)     = $(sech(x))")
    println("csch(x)     = $(csch(x))")
    println("coth(x)     = $(coth(x))")
    println("asinh(x)    = $(asinh(x))")
    println("acosh(x)    = $(acosh(x))")
    println("atanh(x)    = $(atanh(x))")
    println("asech(x)    = $(asech(x))")
    println("acsch(x)    = $(acsch(x))")
    println("acoth(x)    = $(acoth(x))")
    println("sinc(x)     = $(sinc(x))")
    println("deg2rad(x)  = $(deg2rad(x))")
    println("mod2pi(x)   = $(mod2pi(x))")
    println("rad2deg(x)  = $(rad2deg(x))")
    println("abs(x)      = $(abs(x))")
    println("")

    @show x y
    println("")
    println("x + y        = $(x + y)")
    println("x - y        = $(x - y)")
    println("x * y        = $(x * y)")
    println("x / y        = $(x / y)")
    println("x \\ y        = $(x \ y)")
    println("x ^ y        = $(x ^ y)")
    println("atan(x, y)   = $(atan(x, y))")
    println("hypot(x, y)  = $(hypot(x, y))")
    println("log(b, x)    = $(log(b, x)); b = $b")
    println("ldexp(x, n)  = $(ldexp(x, n)); n = $n")
    println("mod(x, y)    = $(mod(x, y))")
    println("rem(x, y)    = $(rem(x, y))")
    println("rem2pi(x, r) = $(rem2pi(x, r)); r = $r")
    # max/min functions are a problem for nansafe jacobian
    # max(0, x) and min(x, 0) seem fine to use
    println("max(x, y)    = $(max(x, y))")
    println("min(x, y)    = $(min(x, y))")
    println("")
    # testing
    println("max_nan(x, y) = $(max_nan(x, y))")
    println("min_nan(x, y) = $(min_nan(x, y))")

    #=
    x = Dual(NaN, 1.0)
    y = Dual(NaN, 0.0)
    z = Dual(NaN, 10.0)
    v = Dual(NaN, 5.0)

    @show max_nan(x,y) max_nan(x,y,z) max_nan(x,y,v,z)
    @show min_nan(x,y) min_nan(x,y,z) min_nan(x,y,v,z)

    @time max_nan(x, y)
    @time max_nan(x, y, z)
    @time max_nan(x, y, z, v)

    @time min_nan(x, y)
    @time min_nan(x, y, z)
    @time min_nan(x, y, z, v)
    =#

    return nothing
end