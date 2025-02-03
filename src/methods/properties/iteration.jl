
abstract type Iteration end
abstract type Implicit <: Iteration end

struct Explicit <: Iteration end        # Runge-Kutta/multistep method is explicit
struct DiagonalImplicit <: Implicit end # Runge-Kutta method is diagonal-implicit
struct FullImplicit <: Implicit end     # Runge-Kutta method is full-implicit
struct SingleImplicit <: Implicit end   # multistep method is implicit
