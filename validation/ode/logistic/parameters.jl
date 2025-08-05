
options = Dict(
    :method => Verner6(),
    # :method => BackwardEuler1(),
    # :method => TrapezoidRuleBDF2(),  # 400.32 k allocations: 35.115 MiB w/ Fixed()
    # :method => AdamsBashforth(; order = 2),
    # :method => AdamsMoulton(; order = 2),
    # :method => Heun2(),

    # :adaptive => Fixed(),
    :adaptive => Embedded(; epsilon = 1e-6, alpha = 1e-6, delta = 1e-6, p_norm = 2.0,
                            pid = PIControl(), limiter = SmoothLimiter(),),
    # :adaptive => Doubling(; epsilon = 1e-6, alpha = 1e-6, delta = 1e-6, p_norm = 2.0,
    #                         pid = PIControl(), limiter = SmoothLimiter(),),
    # :adaptive => CentralDiff(; epsilon = 1e-6, alpha = 1e-6, delta = 1e-6,
    #                            p_norm = 2.0, limiter = SmoothLimiter()),

    # :timer => TimeLimit(; wtime_min = 0),

    :state_jacobian => FiniteJacobian(),
    # :state_jacobian => ForwardJacobian(),

    :root_finder => Newton(; linear_method = RFLUFactorization(),),
    # :root_finder => FixedPoint(),

    :eigenmax => NoEigenMax(),
    # :eigenmax => LinearEigenMax(),
    # :eigenmax => KrylovEigenMax(; krylovdim = 1), # TODO: limit krylovdim to ny

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