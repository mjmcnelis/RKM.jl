module RKMDoubleFloatsExt

import RKM: RKM, RKM_root
import DoubleFloats: IEEEFloat, DoubleFloat

function __init__()
    include("$RKM_root/src/tmp/ext/double_float.jl")
end

end