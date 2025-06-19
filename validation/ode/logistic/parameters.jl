
options = Dict(
    :method => RungeKutta4(),
    # :method => BogackiShampine32(),
    # :method => BackwardEuler1(),
    # :method => TrapezoidRuleBDF2(),  # 400.32 k allocations: 35.115 MiB w/ Fixed()
    # :method => AdamsBashforth(; order = 2),
    # :method => AdamsMoulton(; order = 2),
    # :method => BackwardDifferentiationFormula(; order = 2),   # BDF and NDF currently broken
    # :method => NumericalDifferentiationFormula(; order = 2),
    # :method => HeunEuler21(),

    :adaptive => Fixed(),
    # :adaptive => Embedded(; epsilon = 1e-6, alpha = 1e-6, delta = 1e-6, p_norm = 2.0,
    #                         pid = PIControl(), limiter = SmoothLimiter(),),
    # :adaptive => Doubling(; epsilon = 1e-6, alpha = 1e-6, delta = 1e-6, p_norm = 2.0,
    #                         pid = PIControl(), limiter = SmoothLimiter(),),

    # :timer => TimeLimit(; wtime_min = 0),

    :stage_finder => ImplicitStageFinder(;
                        #  state_jacobian = ForwardJacobian(),
                         state_jacobian = FiniteJacobian(),
                        #  linear_method = LUFactorization(),
                         linear_method = RFLUFactorization(),
                        # TODO: limit krylovdim to ny
                        # eigenmax = KrylovEigenMax(; krylovdim = 1),
                     ),

    :root_finder => Newton(),
    # :root_finder => FixedPoint(),

    :sensitivity => NoSensitivity(),
    # :sensitivity => DecoupledDirect(;
    #                     param_jacobian = FiniteJacobian(),
    #                     # param_jacobian = ForwardJacobian(),
    #                     # jacobian_vector = NaiveJacobianVector(),
    #                     jacobian_vector = FiniteJacobianVector(),
    #                     # jacobian_vector = ForwardJacobianVector(),
    #                 ),

    :interpolator => NoInterpolation(),
    # :interpolator => CubicHermite(),
    # :interpolator => ContinuousFormula(),

    :save_solution => true,
    # :save_solution => false,

    :save_time_derivative => false,
    # :save_time_derivative => true,

    :show_progress => false,
    # :show_progress => true,

    :benchmark_subroutines => false,
    # :benchmark_subroutines => true,

    :precision => Float64
    # :precision => Double64    # sensitivity doesn't work for Double64 rn
    # :precision => BigFloat
);