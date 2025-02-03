
# TODO: determine if butcher table fixed (square matrix) or embedded (not square matrix)

abstract type EmbeddedPair end

struct DefaultPair <: EmbeddedPair end
struct EulerPair <: EmbeddedPair end
struct SecondPair <: EmbeddedPair end
