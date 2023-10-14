# TODO: try adding precision to list
# note: may consider taking out t0, tf out of t_range (repurpose TimeRange as interpolator)
#       could move dt0 to adaptive methods

options = Dict(
    :method => RungeKutta4(),
    # :method => TrapezoidRuleBDF2(),  # 400.32 k allocations: 35.115 MiB w/ Fixed()
    # :method => HeunEuler21(),

    :adaptive => Fixed(),
    # :adaptive => Embedded(; epsilon = 1e-6, p_norm = 2.0),
    # :adaptive => Doubling(; epsilon = 1e-6, p_norm = 2.0),

    :t_range      => TimeRange(; t0 = -10.0, tf = 10.0, dt0 = 1e-4),
    :controller   => TimeStepController(;
                         pid = H312Control(),
                        #  pid = PIControl(),
                         limiter = SmoothLimiter(),
                        #  limiter = PiecewiseLimiter(),
                     ),
    :stage_finder => ImplicitStageFinder(;
                         jacobian_method = ForwardJacobian(),
                        #  jacobian_method = FiniteJacobian(),
                         root_method = Newton(),
                        #  root_method = FixedPoint(),
                         epsilon = 1e-8, max_iterations = 10, p_norm = 2.0,
                     ),
);