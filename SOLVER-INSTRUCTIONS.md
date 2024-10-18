## Repository of instances used for paper 

This repository allows one to create a subset of the instances for the paper from scratch. However, it is much easier to just use the created (and scaled) instances that can be found at [here](https://www.oliverhinder.com/huge-linear-programs).

## Instruction for reproducing Gurobi results

To produce the gurobi results for the paper, we ran 

```{sh}
$ gurobi_cl Crossover=0 Method=[method] [mps_file_name]
```

where method=0 for primal Simplex, 1 for dual and 2 for barrier.
We used Gurobi 11.0.1 (https://www.gurobi.com/downloads/).

Documentation for this method can be found at https://www.gurobi.com/documentation/current/refman/grb_command_line_tool.html


## Instruction for reproducing PDLP results

Download the **source files** for or-tools from the bottom the page of https://developers.google.com/optimization/install. Note we used version 9.10 but we don't expect different behaviour from 9.11.

Build it from source, e.g., run the following command from inside the or-tools directory.

```{sh}
$ cmake -S. -Bbuild -DUSE_COINOR=OFF -DUSE_SCIP=OFF -DBUILD_DEPS:BOOL=O -DCMAKE_C_COMPILER=[gcc_path] -DCMAKE_CXX_COMPILER=[gplusplus_path]
cmake --build build
```

Note that PDLP accepts .mps files but not .mps.gz files. Thus
before using PDLP you need to ungzip the mps files. For example, you can do that via

```{sh}
$ gunzip "xxx.mps.gz"
```

Next, navigate to the directory that you want to run PDLP from,
then run the command:

```{sh}
$ [path to or-tools]/build/bin/pdlp_solve --input [mps filename] --sol_file [solution filename] -- solve_log_file [log file name] \
--params "verbosity_level: 4 num_threads: [number of threads] termination_criteria {detailed_optimality_criteria {eps_optimal_primal_residual_absolute: 1.0e-8 eps_optimal_primal_residual_relative: 0.0 eps_optimal_dual_residual_absolute: 1.0e-8 eps_optimal_dual_residual_relative: 0.0 eps_optimal_objective_gap_absolute: 0.0 eps_optimal_objective_gap_relative: 1.0e-2} eps_primal_infeasible: 1.0e-9 eps_dual_infeasible: 1.0e-9 optimality_norm: OPTIMALITY_NORM_L_INF} use_feasibility_polishing: [true/false] handle_some_primal_gradients_on_finite_bounds_as_residuals: false"
```

This produces the results table. The values for the number of threads
is either 16 or 32 depending on the instance as described in the paper.
The speed up against the number of threads figure can be reproduced by running

```{sh}
$ [path to or-tools]/build/bin/pdlp_solve --input [mps filename] --sol_file [solution filename] -- solve_log_file [log file name] \
--params "verbosity_level: 4 num_threads: [number of threads] termination_criteria {detailed_optimality_criteria {eps_optimal_primal_residual_absolute: 1.0e-8 eps_optimal_primal_residual_relative: 0.0 eps_optimal_dual_residual_absolute: 1.0e-8 eps_optimal_dual_residual_relative: 0.0 eps_optimal_objective_gap_absolute: 0.0 eps_optimal_objective_gap_relative: 1.0e-2} eps_primal_infeasible: 1.0e-9 eps_dual_infeasible: 1.0e-9 optimality_norm: OPTIMALITY_NORM_L_INF iteration_limit: 10000} use_feasibility_polishing: [true/false] handle_some_primal_gradients_on_finite_bounds_as_residuals: false"
```

