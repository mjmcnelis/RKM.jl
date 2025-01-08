
abstract type RootMethod end

struct Newton <: RootMethod end         # TODO: add line search later
struct FixedPoint <: RootMethod end