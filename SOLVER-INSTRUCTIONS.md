## Instances used for the paper

This repository allows one to generate from scratch a subset of the instances for the paper.
See README.md in this repository for instructions on how to generate the instances.

To save time, we recommend downloading the already generated (and scaled) instances from [here](https://www.oliverhinder.com/large-scale-LP-problems).
Download these files to a directory and save this directory path in your linux terminal as a local variable, which will be referenced in example commands below.

```sh
$ path_to_instances=[directory where these instances are]
```

Note that PDLP accepts .mps files but not .mps.gz files. Thus before using PDLP, you need to ungzip the mps.gz files. For example, you can do that via

```sh
$ gunzip "$path_to_instances"/*.mps.gz
```

These files are large, so you may prefer to uncompress instances only as you need them.

## Preparing to run PDLP

Our experiments used OR-Tools version 9.10 (we expect 9.11 to have similar behavior).

To fetch the source OR-Tools version 9.10, run
```sh
git clone --branch v9.10 https://github.com/google/or-tools.git
```

Once you've done that, set

```sh
$ path_to_OR_tools=[path to or-tools]
```

We require a patch to `examples/cpp/CMakeLists.txt` so that `pdlp_solve` is included in the build.
Either apply the patch below or manually delete the indicated lines before proceeding to compile.

```patch
diff --git a/examples/cpp/CMakeLists.txt b/examples/cpp/CMakeLists.txt
index 7fe5b30fcd..2cf8e7e97a 100644
--- a/examples/cpp/CMakeLists.txt
+++ b/examples/cpp/CMakeLists.txt
@@ -50,10 +50,8 @@ list(FILTER CXX_SRCS EXCLUDE REGEX ".*/knapsack_2d_sat.cc")
 list(FILTER CXX_SRCS EXCLUDE REGEX ".*/mps_driver.cc") # crash
 list(FILTER CXX_SRCS EXCLUDE REGEX ".*/multi_knapsack_sat.cc") # crash
 list(FILTER CXX_SRCS EXCLUDE REGEX ".*/network_routing_sat.cc")
-list(FILTER CXX_SRCS EXCLUDE REGEX ".*/pdlp_solve.cc")
 list(FILTER CXX_SRCS EXCLUDE REGEX ".*/pdptw.cc")
 list(FILTER CXX_SRCS EXCLUDE REGEX ".*/shift_minimization_sat.cc")
-list(FILTER CXX_SRCS EXCLUDE REGEX ".*/pdlp_solve.cc")
 list(FILTER CXX_SRCS EXCLUDE REGEX ".*/strawberry_fields_with_column_generation.cc") # Too long
 list(FILTER CXX_SRCS EXCLUDE REGEX ".*/vector_bin_packing_solver.cc")
 list(FILTER CXX_SRCS EXCLUDE REGEX ".*/weighted_tardiness_sat.cc")
 ```

Also, make sure to install any prerequisites for compiling as listed at https://developers.google.com/optimization/install/cpp/source_linux.

Next, compile it, for example, by running the following command from inside the or-tools directory.

```{sh}
$ cd $path_to_OR_tools
$ cmake -S. -Bbuild -DUSE_COINOR=OFF -DUSE_SCIP=OFF -DBUILD_DEPS:BOOL=O -DCMAKE_C_COMPILER=[gcc_path] -DCMAKE_CXX_COMPILER=[gplusplus_path]
$ cmake --build build
```

Check that `build/bin/pdlp_solve` is present after compiling.

## Reproducing Table 3

We will now detail how to reproduce the design-match row of Table 3. The other rows can be produced similarly by modifying the instance_name. Note that we increase num_threads from 16 to 32 for PDLP if it needs more than 256GB of RAM (this only occurs for tsp-gaia-100m).
You may choose to use a different number of threads depending on your machine.
We report the timing recorded by PDLP and Gurobi (and printed to standard output). Note that this excludes the time reading the instance.

To produce the 'without polishing' column:

```{sh}
$ path_to_solution_files=[directory for solution files]
$ path_to_log_files=[directory for log files]
$ instance_name="design-match"
$ "$path_to_OR_tools"/build/bin/pdlp_solve \
--input "$path_to_instances"/"$instance_name".mps \
--sol_file "$path_to_solution_files"/no-polish-"$instance_name".sol \
--solve_log_file "$path_to_log_files"/no-polish-"$instance_name".json \
--params "verbosity_level: 4 num_threads: 16 termination_criteria {detailed_optimality_criteria {eps_optimal_primal_residual_absolute: 1.0e-8 eps_optimal_primal_residual_relative: 0.0 eps_optimal_dual_residual_absolute: 1.0e-8 eps_optimal_dual_residual_relative: 0.0 eps_optimal_objective_gap_absolute: 0.0 eps_optimal_objective_gap_relative: 1.0e-2} eps_primal_infeasible: 1.0e-9 eps_dual_infeasible: 1.0e-9 optimality_norm: OPTIMALITY_NORM_L_INF} use_feasibility_polishing: false handle_some_primal_gradients_on_finite_bounds_as_residuals: false"
```

To produce the 'with polishing' column 

```{sh}
$ "$path_to_OR_tools"/build/bin/pdlp_solve \
--input "$path_to_instances"/"$instance_name".mps \
--sol_file "$path_to_solution_files"/polish-"$instance_name".sol \
--solve_log_file "$path_to_log_files"/polish-"$instance_name".json \
--params "verbosity_level: 4 num_threads: 16 termination_criteria {detailed_optimality_criteria {eps_optimal_primal_residual_absolute: 1.0e-8 eps_optimal_primal_residual_relative: 0.0 eps_optimal_dual_residual_absolute: 1.0e-8 eps_optimal_dual_residual_relative: 0.0 eps_optimal_objective_gap_absolute: 0.0 eps_optimal_objective_gap_relative: 1.0e-2} eps_primal_infeasible: 1.0e-9 eps_dual_infeasible: 1.0e-9 optimality_norm: OPTIMALITY_NORM_L_INF} use_feasibility_polishing: true handle_some_primal_gradients_on_finite_bounds_as_residuals: false"
```

To produce the 'Gurobi barrier' column:

```{sh}
gurobi_cl Crossover=0 Method=2 "$path_to_instances"/"$instance_name".mps.gz  \
LogFile=$path_to_log_files/barrier-"$instance_name".log
```

To produce the 'Gurobi primal simplex' column:

```{sh}
gurobi_cl Crossover=0 Method=0 "$path_to_instances"/"$instance_name".mps.gz \
LogFile=$path_to_log_files/primal-"$instance_name".log
```

To produce the 'Gurobi dual simplex' column:

```{sh}
gurobi_cl Crossover=0 Method=1 "$path_to_instances"/"$instance_name".mps.gz \
LogFile=$path_to_log_files/dual-"$instance_name".log
```

Documentation for gurobi_cl can be found in the [Gurobi command line documentation](https://www.gurobi.com/documentation/current/refman/grb_command_line_tool.html)

## Reproducing Table 4

We get the information on the final primal and dual objectives from the end of the PDLP log files which can be used to reproduce this table.

## Reproducing Figure 1

The speed up against the number of threads figure can be reproduced by running

```{sh}
$ num_threads="4"
$ "$path_to_OR_tools"/build/bin/pdlp_solve \
--input "$path_to_instances"/"$instance_name".mps \
--sol_file "$path_to_solution_files"/thread-test-"$instance_name".sol \
--solve_log_file "$path_to_log_files"/thread-test-"$instance_name".json \
--params "verbosity_level: 4 num_threads: $num_threads termination_criteria {detailed_optimality_criteria {eps_optimal_primal_residual_absolute: 1.0e-8 eps_optimal_primal_residual_relative: 0.0 eps_optimal_dual_residual_absolute: 1.0e-8 eps_optimal_dual_residual_relative: 0.0 eps_optimal_objective_gap_absolute: 0.0 eps_optimal_objective_gap_relative: 1.0e-2} eps_primal_infeasible: 1.0e-9 eps_dual_infeasible: 1.0e-9 optimality_norm: OPTIMALITY_NORM_L_INF iteration_limit: 10000} use_feasibility_polishing: false handle_some_primal_gradients_on_finite_bounds_as_residuals: false"
```

The raw timing values that produced this figure were:

| Number of Threads | Number of Iterations | Time (s)    |
|-------------------|----------------------|-----------|
| 1                 | 10000                | 128863.7  |
| 2                 | 10000                | 65681.4   |
| 4                 | 10000                | 40590.9   |
| 8                 | 10000                | 27838.1   |
| 16                | 10000                | 23296.1   |
| 32                | 10000                | 16974.1   |
