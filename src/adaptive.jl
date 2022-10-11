
abstract type AdaptiveStepSize end

struct Fixed <: AdaptiveStepSize end
struct Doubling <: AdaptiveStepSize end 
struct FiniteDiff <: AdaptiveStepSize end

@kwdef struct Embedded <: AdaptiveStepSize
    pair::EmbeddedPair = DefaultPair()
end
