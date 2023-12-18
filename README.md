Repository for generating large-scale LP test instances for PDLP MPC paper.

Example usage:

```shell
$ mkdir problem-instances
$ julia-1.7 generate-multicommodity-flow.jl --output_file \
    problem-instances/multicommodity-flow-small-test-instance.mps.gz \
    --num_commodities 100 \
    --num_warehouses 30 \
    --num_stores 100
```

For small instances you can also solve the problem using HiGHS and generate
plots of the optimal solution:

```shell
$ mkdir problem-instances
$ julia-1.7 generate-multicommodity-flow.jl --output_file \
    problem-instances/multicommodity-flow-tiny-test-instance.mps.gz \
    --optimize_model \
    --folder_for_plots plot_optimal_flows
```


More examples and documentation will be added in the near future.
