module RKMForwardDiffExt

import RKM: RKM, RKM_root, ForwardJacobian # does this work?
import ForwardDiff: jacobian!, JacobianConfig, DEFAULT_CHUNK_THRESHOLD

include("$RKM_root/src/stage_finder/jacobian/ext/forward_diff_jacobian.jl")

end