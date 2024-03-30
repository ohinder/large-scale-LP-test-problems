mkdir -p small-problem-instances

julia --project generate-multicommodity-flow.jl \
    --output_file small-problem-instances/multicommodity-flow-instance.mps.gz \
    --num_commodities 100 \
    --num_warehouses 30 \
    --num_stores 100 \
    --seed 1

julia --project generate-heat-source-location.jl \
    --output_file small-problem-instances/heat-source-instance-easy.mps.gz \
    --ground_truth_file small-problem-instances/temperature_ground_truth-easy.hdf5 \
    --grid_size 50 \
    --num_source_locations 3 \
    --num_possible_source_locations 100 \
    --num_measurement_locations 200 \
    --seed 1 \
    --maximum_relative_measurement_error 0.0 \
    --pde_solve_tolerance 1e-12

julia --project generate-heat-source-location.jl \
    --output_file small-problem-instances/heat-source-instance-hard.mps.gz \
    --ground_truth_file small-problem-instances/temperature_ground_truth-hard.hdf5 \
    --grid_size 50 \
    --num_source_locations 3 \
    --num_possible_source_locations 100 \
    --num_measurement_locations 50 \
    --seed 2 \
    --maximum_relative_measurement_error 0.0 \
    --pde_solve_tolerance 1e-12

# Note: we calculate epsilon = 1/sqrt(num_treatment_samples) = 0.0141
julia --project design-matching-synthetic.jl \
    --output_file small-problem-instances/synthetic-design-match.mps.gz \
    --epsilon 0.0141 \
    --num_treatment_samples 5000 \
    --num_control_samples 20000 \
    --num_covariates 8 \
    --num_edges_per_treatment 10 \
    --control_shift_magnitude 0.1 \
    --seed 1

julia --project generate-production-inventory.jl \
    --output_file small-problem-instances/production-inventory.mps.gz \
    --num_factories 10 \
    --num_stages 20 \
    --uncertainty_level 0.2
