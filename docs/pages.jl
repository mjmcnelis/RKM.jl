pages = [
    "Home" => "index.md",
    "Overview" => "overview.md",
    "Module" => "module.md",
    "Methods" => ["Runge-Kutta" => ["methods/runge_kutta/runge_kutta.md",
                                    "methods/runge_kutta/tableau.md",
                                    "methods/runge_kutta/debug_table.md",
                                    "Explicit methods" => ["Low-order" => "methods/runge_kutta/explicit/low_order.md",
                                                           "Medium-order" => "methods/runge_kutta/explicit/medium_order.md",
                                                           "High-order" => "methods/runge_kutta/explicit/high_order.md",
                                                           "Very high-order" => "methods/runge_kutta/explicit/very_high_order.md"],
                                    "Implicit methods" => ["Low-order" => "methods/runge_kutta/implicit/low_order.md",
                                                          ]
                                   ],
                 ],
    "Adaptive Time Step" => ["adaptive/constructor.md",
                             "Algorithms" => ["adaptive/algorithms/doubling.md",
                                              "adaptive/algorithms/embedded.md",
                                              "adaptive/algorithms/finite_diff.md"
                                             ],
                            #  "Stepsize control" => ["adaptive/basic_control.md"]
                            ],
    "Solver options" => "solver_options.md",
    "Monitoring" => "monitor/monitor.md",
    "Post-process solution" => ["Solution data" => "solution/solution_data.md",
                                "Dense output" => "solution/dense_output.md",
                                ],
    "Solver statistics" => "statistics.md"
]
