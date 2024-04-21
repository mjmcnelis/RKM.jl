# TODO: try adding precision to list

options = Dict(
    :method => RungeKutta4(),
    # :method => BackwardEuler1(),
    # :method => BackwardDifferentiationFormula2(),
    # :method => TrapezoidRuleBDF2(),  # 400.32 k allocations: 35.115 MiB w/ Fixed()
    # :method => HeunEuler21(),

    :adaptive => Fixed(),
    # :adaptive => Embedded(; epsilon = 1e-6, alpha = 1e-6, delta = 1e-6, p_norm = 2.0),
    # :adaptive => Doubling(; epsilon = 1e-6, alpha = 1e-6, delta = 1e-6, p_norm = 2.0),

    # :timer => TimeLimit(; wtime_min = 0),

    :controller   => TimeStepController(;
                         pid = PIControl(),
                        #  pid = H312Control(),
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
    :save_solution => true,
    # :save_solution => false,

    :show_progress => false
    # :show_progress => true
);