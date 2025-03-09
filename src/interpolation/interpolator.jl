
abstract type Interpolator end
abstract type DenseInterpolator <: Interpolator end

struct NoInterpolation <: Interpolator end
struct CubicHermite <: DenseInterpolator end
struct ContinuousFormula <: DenseInterpolator end
