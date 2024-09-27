module RKMProgressMeterExt

import RKM: RKM, RKM_root
import RKM: create_progress, monitor_progress
import ProgressMeter: Progress, next!

include("$RKM_root/src/progress/ext/progress.jl")

end