
options = Dict(
    # :method => RungeKutta4(),
    # :method => BogackiShampine32(),
    :method => BackwardEuler1(),
    # :method => TrapezoidRuleBDF2(),  # 400.32 k allocations: 35.115 MiB w/ Fixed()
    # :method => AdamsBashforth(; order = 2),
    # :method => AdamsMoulton(; order = 2),
    # :method => BackwardDifferentiationFormula(; order = 2),   # BDF and NDF currently broken
    # :method => NumericalDifferentiationFormula(; order = 2),
    # :method => HeunEuler21(),

    :adaptive => Fixed(),
    # :adaptive => Embedded(; epsilon = 1e-6, alpha = 1e-6, delta = 1e-6, p_norm = 2.0),
    # :adaptive => Doubling(; epsilon = 1e-6, alpha = 1e-6, delta = 1e-6, p_norm = 2.0),

    # :timer => TimeLimit(; wtime_min = 0),

    :controller => TimeStepController(;
                       pid = PIControl(),
                      #  pid = H312Control(),
                       limiter = SmoothLimiter(),
                      #  limiter = PiecewiseLimiter(),
                   ),

    :stage_finder => ImplicitStageFinder(;
                        #  jacobian_method = ForwardJacobian(),
                         jacobian_method = FiniteJacobian(),
                         root_method = Newton(),
                        #  root_method = FixedPoint(),
                        #  linear_method = LUFactorization(),
                         linear_method = RFLUFactorization(),
                         epsilon = 1e-8, max_iterations = 10, p_norm = 2.0,
                     ),

    :sensitivity_method => NoSensitivity(),
    # :sensitivity_method => DecoupledDirect(;
    #                         #    jacobian_method = FiniteJacobian(),
    #                            jacobian_method = ForwardJacobian(),
    #                        ),

    :interpolator => NoInterpolator(),
    # :interpolator => HermiteInterpolator(; dt_save = 0.1),

    :save_solution => true,
    # :save_solution => false,

    :show_progress => false,
    # :show_progress => true,

    :benchmark_subroutines => false,
    # :benchmark_subroutines => true,

    :precision => Float64
    # :precision => Double64
    # :precision => BigFloat
);